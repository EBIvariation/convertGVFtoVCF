import json
import os

import oracledb
from pypika import Query, Table, Schema, Parameter

from ebi_eva_common_pyutils.config import cfg
from ebi_eva_common_pyutils.logger import logging_config as log_cfg

logger = log_cfg.get_logger(__name__)


class DGVaMetadataRetriever:
    """
    The responsibility of this class is to retrieve the metadata for submission to EVA.
    It will extract the required properties to fill in the EVA JSON submission schema
    (https://github.com/EBIvariation/eva-sub-cli/blob/main/eva_sub_cli/etc/eva_schema.json)
    by querying the DGVa database.
    It will generate the metadata submission JSON file.
    """
    def __init__(self, path_to_config_yaml):
        # coming from the config file
        cfg.load_config_file(path_to_config_yaml)  # cfg is a dictionary
        # db connection setup
        self._connection = None
        self._host = self._get_validated_value(cfg, ("DGVA","key_to_host"), str, default_value=None) # get information from the config dictionary
        self._port = self._get_validated_value(cfg, ("DGVA", "key_to_port"), str, default_value=None)
        self._username = self._get_validated_value(cfg, ("DGVA", "key_to_user"), str, default_value=None)
        self._password = self._get_validated_value(cfg, ("DGVA", "key_to_pw"), str, default_value=None)
        print(f"{self._host} {self._port} {self._username} {self._password}")
        # db parameters
        self._max_retries = 3

    def __enter__(self):
        # enables the "with" context manager
        _ = self.connection
        logger.info("open connection")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # closes the connection and resets the state of connection
        if self._connection:
            self._connection.close()
            self._connection = None
            logger.info("Connection closed safely - SUCCESS")

    @property
    def connection(self):
        attempts = 0
        while attempts < self._max_retries:
            # if there is a healthy connection established, return it
            if self._connection is not None and self._connection.is_healthy():
                logger.info("A healthy connection has already been established.")
                return self._connection
            # if there is an unhealthy connection, close it and reset the state
            if self._connection:
                logger.warning("Unhealthy connection. Attempting to reconnect.")
                try:
                    # close the unhealthy connection
                    self._connection.close()
                    logger.info("Connection closed safely.")
                except:
                    pass
                self._connection = None
            # first connection or reconnection
            attempts += 1
            try:
                logger.info("Connecting to the database.")
                self._connection = oracledb.connect(host=self._host, port=self._port,
                                                    user=self._username, password=self._password)
                return self._connection
            except Exception as e:
                logger.error(f"Failed to connect. {e}")
                if attempts >= self._max_retries:
                    logger.error(f"Max connection retries reached: {self._max_retries}.")
                    raise

    def load_from_db(self, query_to_load, query_place_holder_to_load):
        try:
            # create the iterator to process queries
            with self.connection.cursor() as cur:
                # prepare to fetch data
                ###############################################################################
                # WARNING: do not use f-strings or % formatting here, use placeholders instead
                query = query_to_load
                query_placeholders = query_place_holder_to_load
                ###############################################################################
                # run the sql query
                cur.execute(query, query_placeholders)
                # fetch the data from the cursor
                row = cur.fetchall()
                if row:
                    logger.info(f"Fetching metadata query - SUCCESS - {len(dict)} records found")
                    # sql row object is converted to python dict
                    return dict(row)
                else:
                    logger.info("Fetching metadata query - SUCCESS - 0 records found")
                    return {}
        except Exception as e:
            logger.warning(f"Database error: {e}")
            # rollback failed transaction
            self.connection.rollback()
            return {}
        finally:
            self.connection.close()

    def create_json_file(self, json_file_path, study_accession):
        # determine if project new or pre-registered (most projects will be new)
        is_project_preregistered = self._determine_project_pre_registered()
        if is_project_preregistered:
            project_metadata = self._get_project_pre_registered()
        else:
            project_metadata = self._get_project_new()

        # determine if sample new or pre-registered
        is_sample_preregistered = self._determine_sample_pre_registered()
        if is_sample_preregistered:
            sample_metadata = self._get_sample_pre_registered()
        else:
            sample_metadata = self._get_sample_new()

        json_in_eva_format = {
            "submitterDetails": self._get_submitter_details(study_accession),
            "project": project_metadata,
            "analysis": self._get_analysis(),
            "sample": sample_metadata,
            "files": self._get_files()
        }
        with open(json_file_path, 'w') as f:
            json.dump(json_in_eva_format, f, indent=4)
        logger.info(f"Write JSON file for {study_accession}- SUCCESS: {json_file_path}")

    @staticmethod
    def _get_validated_value(config, key_parts_to_get, expected_type, default_value=None):
        # key_parts_to_get is a tuple of multiple values so unpacking
        value_in_config = config.query(*key_parts_to_get, ret_default=default_value)
        if value_in_config is None:
            raise ValueError(f"Missing key required: {key_parts_to_get}")
        if not isinstance(value_in_config, expected_type):
            try:
                # cast the value to its expected type
                value_in_config = expected_type(value_in_config)
            except (ValueError, TypeError):
                raise TypeError(f"Key '{key_parts_to_get}' must be {expected_type.__name__}")
        return value_in_config


    def _get_submitter_details(self, study_accession):
        # create the schema objects
        db = Schema("DGVA")
        # create the table objects
        sc = Table("STUDY_CONTACT", schema=db).as_("sc")
        ds = Table("DGVA_STUDY", schema=db).as_("ds")
        # programmatically create the queries, using named placeholders for readability
        first_name_query = (
            Query.from_(sc)
            .join(ds)
            .on(sc.STUDY_ACCESION == ds.STUDY_ACCESSION)
            .select(sc.FIRST_NAME)
            .where(ds.S_ACC == Parameter(':study_acc'))
        )
        first_name_query_placeholders = {'study_acc': study_accession}
        # submitter detail SQL queries and query placeholders.
        last_name_query = "SELECT * FROM table WHERE x = ? AND y = ? AND z > ?"
        last_name_query_placeholders = ("value3", "value4", study_accession)
        email_query = "SELECT * FROM table WHERE x = ? AND y = ? AND z > ?"
        email_query_placeholders = ("value5", "value6", study_accession)
        centre_query = "SELECT * FROM table WHERE x = ? AND y = ? AND z > ?"
        centre_query_placeholders = ("value7", "value8", study_accession)
        # load metadata
        # expected values in form of {'A': 101, 'B': 102}
        first_name_rows = self.load_from_db(first_name_query.get_sql(), first_name_query_placeholders)
        last_name_rows = self.load_from_db(last_name_query, last_name_query_placeholders)
        email_rows = self.load_from_db(email_query, email_query_placeholders)
        self.laboratory = "PLACEHOLDER VALUE"
        centre_rows = self.load_from_db(centre_query, centre_query_placeholders)
        # store list of submitters
        keys = first_name_rows.keys() # assume keys are the same
        submitter_details_array = []
        for k in keys:
            submitter_object = {
                "lastName": first_name_rows[k],
                "firstName": last_name_rows[k],
                "email": email_rows[k],
                "laboratory": self.laboratory,
                "centre": centre_rows[k]
            }
            submitter_details_array.append(submitter_object)
        return submitter_details_array

    def _get_project_pre_registered(self):
        # return project_object
        # requires a project accession
        # check project accession meets the regex
        pass

    def _get_project_new(self):
        # return project object
        # requries title, description, taxID, centre
        # check taxID is int
        pass

    def _get_analysis(self):
        # return analysis_array
        pass

    def _get_sample_pre_registered(self):
        # return sample_array
        # requires analysisAlias, sampleinVCF, biosample_accession
        pass

    def _get_sample_new(self):
        # return sample_array
        # requires analysisAlias, sampleinVCF, bioSampleObject
        # bioSampleObject requires = name, taxID, scientific_name, release (hold-date) which can be found in DGVA
        # bioSampleObject requires = collection date, geo loc which can be set to unknown/not collected
        pass

    def _get_files(self):
        # return files_array
        # requires analysisAlias, fileName
        pass

    def _determine_project_pre_registered(self):
        # (most projects will be new)
        if self.project_accession:
            is_project_preregistered = True
        else:
            is_project_preregistered = False
        return is_project_preregistered

    def _determine_sample_pre_registered(self):
        if self.biosample_accession:
            is_sample_preregistered = True
        else:
            is_sample_preregistered = False
        return is_sample_preregistered


# def main():
#     # Set up logging functionality
#     ####################################
#     # set path to log and config_file
#
#     #####################################
#     log_cfg.add_file_handler(log_file)
#     logger.info(f"The log file is {log_file}")
#
#     retrieved_dgva_metadata = DGVaMetadataRetriever(config_file)
#     # using with to manage the connection
#     with retrieved_dgva_metadata:
#         retrieved_dgva_metadata.create_json_file(json_file_path="", study_accession="estd22")
#
#
#
# if __name__ == "__main__":
#     main()
