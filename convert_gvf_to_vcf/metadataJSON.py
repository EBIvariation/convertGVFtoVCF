import json
import os
import re
from collections import namedtuple
import oracledb


from pypika import Query, Table, Schema

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
        self._host = self._get_validated_value(cfg, ("DGVA","host"), str, default_value=None) # get information from the config dictionary
        self._port = self._get_validated_value(cfg, ("DGVA", "port"), str, default_value=None)
        self._username = self._get_validated_value(cfg, ("DGVA", "user"), str, default_value=None)
        self._password = self._get_validated_value(cfg, ("DGVA", "password"), str, default_value=None)
        self._service_name = self._get_validated_value(cfg, ("DGVA", "service_name"), str, default_value=None)
        # db parameters
        self._max_retries = 3
    def __enter__(self):
        # enables the "with" context manager
        _ = self.connection
        logger.info("Opening connection using a context manager - SUCCESS")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # closes the connection and resets the state of connection
        if self._connection:
            self._connection.close()
            self._connection = None
            logger.info("Closing connection safely using a context manager - SUCCESS")

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
                                                    user=self._username, password=self._password, service_name=self._service_name)
                return self._connection
            except Exception as e:
                logger.error(f"Failed to connect. {e}")
                if attempts >= self._max_retries:
                    logger.error(f"Max connection retries reached: {self._max_retries}.")
                    raise

    # def load_from_db(self, query_to_load, query_place_holder_to_load):
    def load_from_db(self, query_to_load):
        try:
            # create the iterator to process queries
            with self.connection.cursor() as cur:
                # prepare to fetch data
                ###############################################################################
                # WARNING: do not use f-strings or % formatting here, use placeholders instead
                query = query_to_load
                # query_placeholders = query_place_holder_to_load
                ###############################################################################
                # run the sql query
                # cur.execute(query, query_placeholders)
                logger.info(f"Executing the following query: {query}")
                cur.execute(query)
                # fetch the data from the cursor
                logger.info(f"Fetching the data....")
                row = cur.fetchall()
                logger.info(f" row fetched: {row}")
                if row:
                    logger.info(f"Fetching metadata query - SUCCESS - {len(row)} records found")
                    # sql row object is converted to python dict
                    row_dict = dict(enumerate(row))
                    # row_dict = [{i: val[0] for i, val in enumerate(row)}]
                    return row_dict
                else:
                    logger.info("Fetching metadata query - SUCCESS - 0 records found")
                    return {}
        except Exception as e:
            logger.warning(f"Database error: {e}")
            logger.warning("Rolling back")
            # rollback failed transaction
            self.connection.rollback()
            return {}
        # finally:
        #     self.connection.close()

    def create_json_file(self, json_file_path, study_accession, vcf_output):
        # determine if project new or pre-registered (most projects will be new)
        is_project_preregistered, project_accession = self._determine_project_pre_registered(study_accession)
        if is_project_preregistered:
            project_metadata = self._get_project_pre_registered(project_accession)
        else:
            project_metadata = self._get_project_new(study_accession)

        # determine if sample new or pre-registered
        # list of tuples [(is_sample_preregistered, biosample_accession)]
        sample_registration_statuses = self._determine_sample_pre_registered(study_accession)
        sample_metadata_array = []
        for sample_status in sample_registration_statuses:
            if sample_status.is_sample_preregistered:
                sample_metadata = self._get_sample_pre_registered(study_accession, sample_status.accession)
                sample_metadata_array.append(sample_metadata)
            else:
                sample_metadata = self._get_sample_new(study_accession)
                sample_metadata_array.append(sample_metadata)
        json_in_eva_format = {
            "submitterDetails": self._get_submitter_details(study_accession),
            "project": project_metadata,
            "analysis": self._get_analysis(study_accession),
            "sample": sample_metadata_array,
            "files": self._get_files(study_accession, vcf_output)
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

    # THESE GETTERS GET THE RELEVANT JSON OBJECT
    def _get_submitter_details(self, study_accession):
        """ Formats the submitter details from the DGVa database. Creates a submitter object for each submitter.
        :param: study_accession - expected format ^(estd|nstd)\d+$
        :return: submitter_details_array
        """
        logger.info("Fetching submitter details.")
        submitter_details_query = self._fetch_submitter_details(study_accession)
        # load metadata
        # expected values in form of DICTIONARY OF ENUMERATED TUPLES {INDEX: (TUPLE)}
        # {0: ("FIRST_NAME", "LAST_NAME", "EMAIL", "CENTRE"), 1: ("FIRST_NAME", "LAST_NAME", "EMAIL", "CENTRE")}
        submitter_details_dict = self.load_from_db(submitter_details_query.get_sql(quote_char=None))
        # store list of submitters
        keys = submitter_details_dict.keys() # assume keys are the same
        submitter_details_array = []
        for key in keys:
            submitter_object = {
                "lastName": submitter_details_dict[key][1] or "UNSPECIFIED-LASTNAME",
                "firstName": submitter_details_dict[key][0] or "UNSPECIFIED-FIRSTNAME",
                "email": submitter_details_dict[key][2] or "UNSPECIFIED-EMAIL",
                "laboratory": "UNSPECIFIED-LABORATORY", # because this is not found in dgva
                "centre": submitter_details_dict[key][3] or "UNSPECIFIED-CENTRE"
            }
            for eva_field_name,value in submitter_object.items():
                if "UNSPECIFIED" in str(value):
                    logger.info(f"Fetching {eva_field_name} - FAILURE - {eva_field_name} not found: {value}.")
                else:
                    logger.info(f"Fetching {eva_field_name} - SUCCESS - {eva_field_name} found.")
            submitter_details_array.append(submitter_object)

        return submitter_details_array

    def _get_project_pre_registered(self, project_accession):
        # check project accession meets the regex
        project_accession_pattern = r"^PRJ(E|D|N)[A-Z][0-9]+$"
        assert re.fullmatch(project_accession_pattern, project_accession), f"Project accession does not meet the regex: {project_accession_pattern}"
        project_object = {
            "projectAccesion" : project_accession
        }
        return project_object

    def _get_project_new(self, study_accession):
        MAX_PROJECT_DESCRIPTION_LENGTH = 5000
        MAX_PROJECT_TITLE_LENGTH = 500
        # programmatically create the queries
        # PROJECT TITLE
        project_title = self._fetch_project_title(study_accession)
        # PROJECT DESCRIPTION
        project_description = self._fetch_project_description(study_accession)
        # PROJECT TAX ID
        project_tax_id = self._fetch_tax_id(study_accession)
        #PROJECT CENTRE
        project_centre = self._fetch_centre(study_accession)
        # performing checks
        assert len(project_description) <= MAX_PROJECT_DESCRIPTION_LENGTH, f"Project description exceeded length: {MAX_PROJECT_DESCRIPTION_LENGTH}"
        assert len(project_title) <= MAX_PROJECT_TITLE_LENGTH, f"Project title exceeded length: {MAX_PROJECT_TITLE_LENGTH}"
        assert isinstance(project_tax_id, int), f"Project Tax ID must be an int: {project_tax_id} is {type(project_tax_id)}"

        #TODO: add non-EVA-required metadata to improve metadata completeness.

        # not required: publications, parentProject, childProject, peerProjects, hold-date, links
        # check publications meet regex"pattern": "^[^:,]+?:[^:,]+?$"
        # check parent/child/peer projects meet regex

        # required: title, description, taxID, centre
        project_object = {
            "title": project_title,
            "description": project_description,
            "taxID": project_tax_id,
            "centre": project_centre
        }
        return project_object


    def _get_analysis(self, study_accession):
        # return analysis_array
        # required: analysisTitle, analysisAlias, description, experimentType, reference_genome
        logger.info("Fetching Analysis details.")
        analysis_analysis_alias = self._fetch_analysis_alias(study_accession)
        analysis_analysis_title = "UNSPECIFIED-TITLE"
        analysis_analysis_description = self._fetch_analysis_description(study_accession)
        analysis_experiment_type = self._fetch_experiment_type(study_accession)
        # TODO: DISCUSSION: turning reference into genbank accessions: the reference in DGVA is not GCA (example values: NCBI35, MGSCv37). Also, how to handle null value here?
        analysis_reference_genome = self._fetch_reference_genome(study_accession)
        analysis_array = []
        analysis_object = {
            "analysisTitle": analysis_analysis_title,
            "analysisAlias": analysis_analysis_alias,
            "description": analysis_analysis_description,
            "experimentType": analysis_experiment_type,
            "referenceGenome": analysis_reference_genome
        }
        placeholder_keys = [k for k,v in analysis_object.items() if v is None or "UNSPECIFIED" in v]
        for placeholder_key in placeholder_keys:
            logger.info(f"{placeholder_key} not found. Adding placeholder: {analysis_object[placeholder_key]}")
        analysis_array.append(analysis_object)
        return analysis_array

    def _get_sample_pre_registered(self, study_accession, biosample_accession):
        # requires analysisAlias, sampleinVCF, biosample_accession
        sample_analysis_alias = self._fetch_analysis_alias(study_accession)
        sample_sampleinvcf = "UNSPECIFIED_SAMPLE_IN_VCF"
        sample_object = {
            "analysisAlias": [].append(sample_analysis_alias),
            "sampleInVCF": sample_sampleinvcf,
            "bioSampleAccession": biosample_accession
        }
        return sample_object

    def _get_sample_new(self, study_accession):
        # return sample object
        # requires analysisAlias, sampleinVCF, bioSampleObject
        # bioSampleObject requires = name, taxID, scientific_name, release (hold-date) which can be found in DGVA
        # bioSampleObject requires = collection date, geo loc which can be set to unknown/not collected
        sample_analysis_alias = self._fetch_analysis_alias(study_accession)
        sample_sampleinvcf = "UNSPECIFIED_SAMPLE_IN_VCF"
        # TODO: fix the error, as there multiple sample names,
        sample_name = self._fetch_sample_id(study_accession)
        sample_tax_id = self._fetch_tax_id(study_accession)
        scientific_name = self._fetch_scientific_name(study_accession)
        hold_date = self._fetch_hold_date(study_accession)
        collection_date = "not provided"
        geographic_location_country_and_or_sea = "not provided"
        biosample_object = {
            "sample_title": sample_name,
            "scientific_name": scientific_name,
            "tax_id": sample_tax_id,
            "collection date": collection_date,
            "geographic location (country and/or sea)": geographic_location_country_and_or_sea
        }
        sample_object = {
            "analysisAlias": sample_analysis_alias,
            "sampleInVCF": sample_sampleinvcf,
            "bioSampleObject": biosample_object
        }
        return sample_object

    def _get_files(self, study_accession, vcf_output):
        # return files_array
        # requires analysisAlias, fileName
        files_analysis_alias = self._fetch_analysis_alias(study_accession)
        files_file_name = self._fetch_file_name(vcf_output)
        files_array = []
        files_object = {
            "analysisAlias": files_analysis_alias,
            "fileName": files_file_name
        }
        files_array.append(files_object)
        return files_array

    # THESE DETERMINE IF NEW OR PRE_REG
    def _determine_project_pre_registered(self, study_accession):
        """ Determines if the project is new or pre-registered.
        :param: study_accession - expected format ^(estd|nstd)\d+$
        :return is_project_preregistered, project_accession : boolean and string format ^(PRJ)[A-Z]{2}\d+$
        """
        # create the schema objects
        db = Schema("DGVA")
        # create the table objects
        ds = Table("DGVA_STUDY", schema=db).as_("ds")
        # create the query
        project_accession_query = (Query
                                   .from_(ds)
                                   .select(ds.BIOPROJECT_ACCESSION)
                                   .where(ds.STUDY_ACCESSION == study_accession)
                                   )
        project_accession_dict = self.load_from_db(project_accession_query.get_sql(quote_char=None))
        is_project_preregistered = False
        for value in project_accession_dict.values():
            if isinstance(value, tuple):
                if not None in value:
                    is_project_preregistered = True
                    project_accession = next(iter(project_accession_dict.values()))[0]
                    logger.info(f"Determining if the project is pre-registered - SUCCESS - Project found: {project_accession}.")
                    return is_project_preregistered, project_accession
                else:
                    logger.info(f"Determining if the project is pre-registered - FAILURE - No project found.")
                    return is_project_preregistered, None

    def _determine_sample_pre_registered(self, study_accession):
        # create the schema objects
        db = Schema("DGVA")
        # create the table objects
        dsamp = Table("DGVA_SAMPLE", schema=db).as_("dsamp")
        # create the query
        sample_accession_query = (Query
                                   .from_(dsamp)
                                   .select(dsamp.BIOSAMPLE_ACCESSION)
                                   .where(dsamp.STUDY_ACCESSION == study_accession)
                                   )
        sample_accession_dict = self.load_from_db(sample_accession_query.get_sql(quote_char=None))
        sample_accession_and_status_list = []
        SampleStatus = namedtuple('SampleStatus', ['is_sample_preregistered', 'sample_accession'])
        is_sample_preregistered = False
        for value in sample_accession_dict.values():
            if isinstance(value, tuple):
                if not None in value:
                    # DGVa does not have any samples with a BioSample accession. We do not expect this to be populated.
                    is_sample_preregistered = True
                    sample_accession = next(iter(sample_accession_dict.values()))[0]
                    sample_accession_and_status_list.append(SampleStatus(is_sample_preregistered, sample_accession))
                    # return is_sample_preregistered, sample_accession
                    logger.info(f"Determining if sample is pre-registered - SUCCESS - Sample found: {SampleStatus(is_sample_preregistered, sample_accession)}.")
                else:
                    sample_accession_and_status_list.append(SampleStatus(is_sample_preregistered, None))
                    # return is_sample_preregistered, None
                    logger.info(f"Determining if sample is pre-registered - FAILURE - Sample not found.")
        return sample_accession_and_status_list

    # VALIDATING FETCH RESULTS OR USING PLACEHOLDER
    def validate_fetch_result(self, eva_field_name, fetch_result_dict):
        if fetch_result_dict != {}:
            fetch_result = next(iter(fetch_result_dict.values()))[0]
            logger.info(f"Fetching {eva_field_name} - SUCCESS - {eva_field_name} found: {fetch_result}.")
        else:
            fetch_result = f"UNSPECIFIED_{eva_field_name}"
            logger.info(f"Fetching {eva_field_name}  - FAILURE - {eva_field_name} not found. Adding placeholder.")
        return fetch_result

    #### THESE FETCH FIELDS FROM THE DB
    # SUBMITTER DETAILS SECTION
    def _fetch_submitter_details(self, study_accession):
        # create the schema objects
        db = Schema("DGVA")
        # create the table objects
        sc = Table("STUDY_CONTACT", schema=db).as_("sc")
        ds = Table("DGVA_STUDY", schema=db).as_("ds")
        # programmatically create the queries
        submitter_details_query = (
            Query.from_(sc)
            .join(ds)
            .on(sc.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .select(sc.FIRST_NAME, sc.LAST_NAME, sc.CONTACT_EMAIL, sc.AFFILIATION_NAME)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        return submitter_details_query
    # PROJECT SECTION
    def _fetch_project_title(self, study_accession):
        # create the schema objects
        db = Schema("DGVA")
        # create the table objects
        ds = Table("DGVA_STUDY", schema=db).as_("ds")
        sub = Table("DGVA_SUBMISSION", schema=db).as_("sub")
        project_title_query = (
            Query.from_(ds)
            .join(sub)
            .on(ds.STUDY_ACCESSION == sub.STUDY_ACCESSION)
            .select(ds.DISPLAY_NAME)
            .where(ds.STUDY_ACCESSION == study_accession)

        )
        project_title_dict = self.load_from_db(project_title_query.get_sql(quote_char=None))
        project_title = next(iter(project_title_dict.values()))[0]
        return project_title

    def _fetch_project_description(self, study_accession):
        # create the schema objects
        db = Schema("DGVA")
        # create the table objects
        stc = Table("STUDY_TYPE_CV", schema=db).as_("stc")
        ds = Table("DGVA_STUDY", schema=db).as_("ds")
        project_description_query = (
            Query.from_(ds)
            .join(stc)
            .on(ds.STUDY_TYPE == stc.STUDY_TYPE)
            .select(ds.STUDY_DESCRIPTION)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        project_description_dict = self.load_from_db(project_description_query.get_sql(quote_char=None))
        project_description = next(iter(project_description_dict.values()))[0]
        return project_description

    def _fetch_tax_id(self, study_accession):
        # create the schema objects
        db = Schema("DGVA")
        # create the table objects
        ds = Table("DGVA_STUDY", schema=db).as_("ds")
        so = Table("STUDY_ORGANISM", schema=db).as_("so")
        tax_id_query = (
            Query.from_(ds)
            .join(so)
            .on(so.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .select(so.TAXONOMY_ID)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        tax_id_dict = self.load_from_db(tax_id_query.get_sql(quote_char=None))
        tax_id = next(iter(tax_id_dict.values()))[0]
        return tax_id

    def _fetch_centre(self, study_accession):
        # create the schema objects
        db = Schema("DGVA")
        # create the table objects
        ds = Table("DGVA_STUDY", schema=db).as_("ds")
        sc = Table("STUDY_CONTACT", schema=db).as_("sc")
        project_centre_query = (
            Query.from_(sc)
            .join(ds)
            .on(sc.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .select(sc.AFFILIATION_NAME)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        project_centre_dict = self.load_from_db(project_centre_query.get_sql(quote_char=None))
        project_centre = next(iter(project_centre_dict.values()))[0]
        return project_centre
    # ANALYSIS SECTION
    def _fetch_analysis_alias(self, study_accession):
        # create the schema objects
        db = Schema("DGVA")
        # create the table objects
        ds = Table("DGVA_STUDY", schema=db).as_("ds")
        sa = Table("STUDY_ALIAS", schema=db).as_("sa")
        # ANALYSIS ALIAS
        analysis_alias_query = (
            Query.from_(ds)
            .join(sa)
            .on(sa.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .select(sa.STUDY_ALIAS_NAME)
            .where(ds.STUDY_ACCESSION == study_accession)

        )
        analysis_alias_dict = self.load_from_db(analysis_alias_query.get_sql(quote_char=None))
        analysis_alias = self.validate_fetch_result("analysisAlias", analysis_alias_dict)
        return analysis_alias

    def _fetch_analysis_description(self, study_accession):
        # create the schema objects
        db = Schema("DGVA")
        # create the table objects
        ds = Table("DGVA_STUDY", schema=db).as_("ds")
        dss = Table("DGVA_SAMPLESET", schema=db).as_("dss")
        da = Table("DGVA_ANALYSIS", schema=db).as_("da")
        rsa = Table("REFERENCE_SAMPLESET_ANALYSIS", schema=db).as_("rsa")
        analysis_description_query = (
            Query.from_(dss)
            .join(ds).on(dss.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .join(rsa).on(dss.SAMPLESET_ID == rsa.SAMPLESET_ID)
            .join(da).on(rsa.ANALYSIS_ID == da.ANALYSIS_ID)
            .select(da.ANALYSIS_DESCRIPTION)
            .where(dss.STUDY_ACCESSION == study_accession)
        )
        analysis_description_dict = self.load_from_db(analysis_description_query.get_sql(quote_char=None))
        analysis_description = self.validate_fetch_result("description", analysis_description_dict)
        return analysis_description

    def _fetch_experiment_type(self, study_accession):
        # create the schema objects
        db = Schema("DGVA")
        # create the table objects
        de = Table("DGVA_EXPERIMENT", schema=db).as_("de")
        ds = Table("DGVA_STUDY", schema=db).as_("ds")
        experiment_type_query = (
            Query.from_(de)
            .join(ds).on(de.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .select(de.EXPERIMENT_TYPE)
            .where(de.STUDY_ACCESSION == study_accession)
        )
        experiment_type_dict = self.load_from_db(experiment_type_query.get_sql(quote_char=None))
        experiment_type = next(iter(experiment_type_dict.values()))[0]
        return experiment_type

    def _fetch_reference_genome(self, study_accession):
        # TODO: reference genome is not in GCA needs converting
        pass
    # SAMPLE SECTION
    def _fetch_sample_id(self, study_accession):
        # create the schema objects
        db = Schema("DGVA")
        # create the table objects
        dsamp = Table("DGVA_SAMPLE", schema=db).as_("dsamp")
        sample_id_query = (
            Query.from_(dsamp)
            .select(dsamp.SUBMITTER_SAMPLE_ID)
            .where(dsamp.STUDY_ACCESSION == study_accession)
        )
        sample_id_dict = self.load_from_db(sample_id_query.get_sql(quote_char=None))
        sample_id = next(iter(sample_id_dict.values()))[0]
        return sample_id

    def _fetch_hold_date(self, study_accession):
        # create the schema objects
        db = Schema("DGVA")
        # create the table objects
        ds = Table("DGVA_STUDY", schema=db).as_("ds")
        hold_date_query = (
            Query.from_(ds)
            .select(ds.HOLD_DATE)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        hold_date_dict = self.load_from_db(hold_date_query.get_sql(quote_char=None))
        hold_date = next(iter(hold_date_dict.values()))[0]
        return hold_date

    def _fetch_scientific_name(self, study_accession):
        # create the schema objects
        db = Schema("DGVA")
        # create the table objects
        so = Table("STUDY_ORGANISM", schema=db).as_("so")
        ds = Table("DGVA_STUDY", schema=db).as_("ds")
        oc = Table("ORGANISM_CV", schema=db).as_("oc")
        scientific_name_query = (
            Query.from_(ds)
            .join(so).on(so.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .join(oc).on(so.TAXONOMY_ID == oc.TAXONOMY_ID)
            .select(oc.SPECIES_LATIN_NAME)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        scientific_name_dict = self.load_from_db(scientific_name_query.get_sql(quote_char=None))
        scientific_name = next(iter(scientific_name_dict.values()))[0]
        return scientific_name
    # FILES SECTION
    def _fetch_file_name(self, vcf_output):
        file_name = os.path.basename(vcf_output)
        return file_name
