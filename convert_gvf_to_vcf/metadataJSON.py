import json
import re

import psycopg2
import yaml

from ebi_eva_common_pyutils.config import cfg
from ebi_eva_common_pyutils.logger import logging_config as log_cfg

logger = log_cfg.get_logger(__name__)


class MetadataManager:
    """
    The responsibility of this class is to manage the metadata for submission to EVA.
    It will extract the required properties to fill in the EVA JSON submission schema
    (https://github.com/EBIvariation/eva-sub-cli/blob/main/eva_sub_cli/etc/eva_schema.json)
    by querying the DGVa database.
    It will generate the metadata submission JSON file.
    """
    def __init__(self, path_to_config_yaml, path_to_sql_queries_yaml):
        # coming from the config file
        cfg.load_config_file(path_to_config_yaml)  # cfg is a dictionary
        # db connection setup
        self.connection = None # no connection yet
        self.host = cfg.get("key_to_host") # get information from the config dictionary
        self.port = cfg.get("key_to_port")

        # load sql queries and their placeholders
        sql_map = self.load_sql_registry(path_to_sql_queries_yaml)
        # coming from DGVa
        # the following are dict values obtained from the DGVa database, may need to parse to obtain the value itself
        # SUBMITTER DETAILS
        self.first_name = self.load_from_db(sql_map["SUBMITTER_DETAILS"]["first_name_query"], sql_map["SUBMITTER_DETAILS"]["first_name_query_placeholders"])
        self.last_name = self.load_from_db(sql_map["SUBMITTER_DETAILS"]["last_name_query"], sql_map["SUBMITTER_DETAILS"]["last_name_query_placeholders"])
        self.email = self.load_from_db(sql_map["SUBMITTER_DETAILS"]["email_query"], sql_map["SUBMITTER_DETAILS"]["email_query_placeholders"])
        self.laboratory = "PLACEHOLDER VALUE"
        self.centre = self.load_from_db(sql_map["SUBMITTER_DETAILS"]["centre_query"], sql_map["SUBMITTER_DETAILS"]["centre_query_placeholders"])
        # PROJECT - PREREGISTERED
        self.project_accession = self.load_from_db(sql_map["PROJECT_PREREGISTERED"]["project_accession_query"],
                                                   sql_map["PROJECT_PREREGISTERED"]["project_accession_query_placeholders"]) # required
        regex_pattern = r"^PRJ(E|D|N)[A-Z][0-9]+$" # regex for project accession
        assert re.fullmatch(regex_pattern, self.project_accession), f"Invalid project accession: {self.project_accession} "
        # PROJECT - NEW
        self.title = self.load_from_db(sql_map["PROJECT_NEW"]["title_query"], sql_map["PROJECT_NEW"]["title_query_placeholders"])
        self.project_description = self.load_from_db(sql_map["PROJECT_NEW"]["project_description_query"], sql_map["PROJECT_NEW"]["project_description_query_placeholders"])
        self.tax_id = self.load_from_db(sql_map["PROJECT_NEW"]["tax_id_query"],
            sql_map["PROJECT_NEW"]["tax_id_query_placeholders"])
        assert isinstance(self.tax_id, int), f"Taxonomy id must be an integer. Taxa ID provided: {self.tax_id}"
        # ANALYSIS
        # SAMPLES - PREREGISTERED
        # SAMPLES - NEW
        # FILES

    @staticmethod
    def load_sql_registry(path_to_sql_queries_yaml):
        with open(path_to_sql_queries_yaml, "r") as query_file:
            return yaml.safe_load(query_file)

    def load_from_db(self, query_to_load, query_place_holder_to_load):
        try:
            # create a connection to the database
            self.connection = psycopg2.connect(host=self.host, port=self.port, database=self.database,
                                               username=self.username, password=self.pw)
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
            logger.info("Connection closed safely - SUCCESS")

    def get_submitter_details(self):
        submitter_details_array = [
            {
                "lastName": self.last_name,
                "firstName": self.first_name,
                "email": self.email,
                "laboratory": self.laboratory,
                "centre": self.centre
            }
        ]
        return submitter_details_array

    def get_project_pre_registered(self):
        # return project_object
        # requires a project accession
        # check project accession meets the regex
        pass

    def get_project_new(self):
        # return project object
        # requries title, description, taxID, centre
        # check taxID is int
        pass

    def get_analysis(self):
        # return analysis_array
        pass

    def get_sample_pre_registered(self):
        # return sample_array
        # requires analysisAlias, sampleinVCF, biosample_accession
        pass

    def get_sample_new(self):
        # return sample_array
        # requires analysisAlias, sampleinVCF, bioSampleObject
        # bioSampleObject requires = name, taxID, scientific_name, release (hold-date) which can be found in DGVA
        # bioSampleObject requires = collection date, geo loc which can be set to unknown/not collected
        pass

    def get_files(self):
        # return files_array
        # requires analysisAlias, fileName
        pass


    def determine_project_pre_registered(self):
        if self.project_accession:
            is_project_preregistered = True
        else:
            is_project_preregistered = False
        return is_project_preregistered

    def determine_sample_pre_registered(self):
        if self.biosample_accession:
            is_sample_preregistered = True
        else:
            is_sample_preregistered = False
        return is_sample_preregistered

    def create_json_file(self):
        # determine if project new or pre-registered
        is_project_preregistered = self.determine_project_pre_registered()
        if is_project_preregistered:
            project_metadata = self.get_project_pre_registered()
        else:
            project_metadata = self.get_project_new()

        # determine if sample new or pre-registered
        is_sample_preregistered = self.determine_sample_pre_registered()
        if is_sample_preregistered:
            sample_metadata = self.get_sample_pre_registered()
        else:
            sample_metadata = self.get_sample_new()

        json_in_eva_format = {
            "submitterDetails": self.get_submitter_details(),
            "project": project_metadata,
            "analysis": self.get_analysis(),
            "sample": sample_metadata,
            "files": self.get_files()
        }
        return json_in_eva_format

    def write_json_file(self, json_file_path, json_in_eva_format):
        # provide it with the dictionary and write the file using json.dump
        with open(json_file_path, 'w') as f:
            json.dump(json_in_eva_format, f, indent=4)
        logger.info(f"Write JSON file - SUCCESS: {json_file_path}")

# def main():
    #manager = MetadataManager(path_to_query_mapper, path_to_sql_queries_yaml)
    # root_json = manager.create_json_file()
    # manager.write_json_file(path_to_json_file, json_in_eva_format)

#
# if __name__ == "__main__":
#     main()
