import hashlib
import json
import os
import re
from datetime import datetime
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
                row = cur.fetchall()
                logger.info(f"Fetching the data: {row}")
                if row:
                    logger.info(f"Fetching metadata query - SUCCESS - {len(row)} records found")
                    # sql row object is a list of tuples converted to python dict
                    row_dict = dict(enumerate(row))
                    # row_dict = [{i: val[0] for i, val in enumerate(row)}]
                    return row_dict
                else:
                    logger.info("Fetching metadata query - SUCCESS - 0 records found")
                    return {}
        except Exception as e:
            logger.warning(f"Database error: {e}")
            return {}

    def create_json_file(self, json_file_path, study_accession, vcf_output, assembly, assembly_report):
        # determine if project new or pre-registered (most projects will be new)
        is_project_preregistered, project_accession = self._determine_project_pre_registered(study_accession)
        if is_project_preregistered:
            project_metadata = self._get_project_pre_registered(project_accession)
        else:
            project_metadata = self._get_project_new(study_accession)

        # determine if sample new or pre-registered
        # list of tuples [(is_sample_preregistered, biosample_accession)]
        sample_ids = self._fetch_sample_id_list(study_accession)
        sample_registration_statuses = self._determine_sample_pre_registered(study_accession)

        sample_metadata_array = []

        for sample_status in sample_registration_statuses:
            if sample_status.is_sample_preregistered:
                sample_metadata = self._get_sample_pre_registered(study_accession, sample_status.sample_accession, sample_status.sample_id)
                sample_metadata_array.append(sample_metadata)
            else:
                # assuming all samples are not pre-registered
                sample_metadata = self._get_sample_new(study_accession, sample_status.sample_id)
                sample_metadata_array.append(sample_metadata)
        json_in_eva_format = {
            "submitterDetails": self._get_submitter_details(study_accession),
            "project": project_metadata,
            "analysis": self._get_analysis(study_accession, vcf_output, assembly, assembly_report),
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

    # THESE GETTERS GET THE RELEVANT JSON OBJECT FOR THE SECTION
    def _get_submitter_details(self, study_accession):
        """ Formats the submitter details from the DGVa database. Creates a submitter object for each submitter.
        :param: study_accession - expected format ^(estd|nstd)\d+$
        :return: submitter_details_array
        """
        logger.info("Fetching submitter details.")
        submitter_details_dict = self._fetch_submitter_details_dictionary(study_accession)
        # load metadata
        # expected values in form of DICTIONARY OF ENUMERATED TUPLES {INDEX: (TUPLE)}
        # {0: ("FIRST_NAME", "LAST_NAME", "TELEPHONE", "EMAIL", "CENTRE", "ADDRESS"), 1: ("FIRST_NAME", "LAST_NAME", "TELEPHONE", "EMAIL", "CENTRE", "ADDRESS")}
        # store list of submitters
        keys = submitter_details_dict.keys() # assume keys are the same
        submitter_details_array = []
        for key in keys:
            submitter_object = {
                "lastName": submitter_details_dict[key][1] or "",
                "firstName": submitter_details_dict[key][0] or "",
                "telephone": submitter_details_dict[key][2] or "",
                "email": submitter_details_dict[key][3] or "",
                "laboratory": "", # expected to be empty because this is not found in dgva
                "centre": submitter_details_dict[key][4] or "",
                "address": submitter_details_dict[key][5] or ""
            }
            for eva_field_name,value in submitter_object.items():
                if str(value) == "":
                    logger.error(f"Fetching {eva_field_name} - FAILURE - {eva_field_name} not found: Value set to empty string.")
                else:
                    logger.info(f"Fetching {eva_field_name} - SUCCESS - {eva_field_name} found.")
            submitter_details_array.append(submitter_object)

        return submitter_details_array

    def _get_project_pre_registered(self, project_accession):
        # check project accession meets the regex
        project_accession_pattern = r"^PRJ(E|D|N)[A-Z][0-9]+$"
        assert re.fullmatch(project_accession_pattern, project_accession), f"Project accession {project_accession} does not meet the regex: {project_accession_pattern}"
        logger.info(f"Project accession {project_accession} has been found as pre-registered.")
        project_object = {
            "projectAccession" : project_accession
        }
        return project_object

    def _get_project_new(self, study_accession):
        logger.info("Fetching Project details.")
        project_title = self._fetch_project_title(study_accession)
        project_description = self._fetch_project_description(study_accession)
        project_tax_id = self._fetch_tax_id(study_accession)
        project_centre_per_submitter = self._fetch_centre(study_accession)
        project_centre = project_centre_per_submitter[0] # choose the first submitter's centre

        project_publications = self._fetch_project_publications(study_accession)
        project_parent_project = self._fetch_project_parent_project(study_accession)
        project_child_project = "" # expect no child projects to exist
        project_peer_project = "" # expect no child projects to exist
        project_links = self._fetch_project_links(study_accession) # expect this to be only one link
        project_hold_date = self._fetch_hold_date(study_accession)

        # formatting
        project_hold_date, project_links, project_parent_project, pubmed_publications = self._format_project(
            project_hold_date, project_links, project_parent_project, project_publications)
        # ensure compliance with EVA JSON schema https://github.com/EBIvariation/eva-sub-cli/blob/main/eva_sub_cli/etc/eva_schema.json
        self.validate_project(project_description, project_hold_date, project_parent_project,
                              project_tax_id, project_title, pubmed_publications)

        logger.info(f"Project accession has not been found. Creating a new project. ")
        # required: title, description, taxID, centre
        # not required: publications, parentProject, childProject, peerProjects, hold-date, links
        project_object_all = {
            "title": project_title,
            "description": project_description,
            "taxId": project_tax_id,
            "centre": project_centre,
            "publications": pubmed_publications,
            "parentProject": project_parent_project,
            "childProject": project_child_project,
            "peerProject": project_peer_project,
            "links": project_links,
            "holdDate": project_hold_date
        }
        project_object = {k:v for k,v in project_object_all.items() if v}
        return project_object



    def _get_analysis(self, study_accession, vcf_output, assembly, assembly_report):
        # return analysis_array
        # required: analysisTitle, analysisAlias, description, experimentType, reference_genome
        logger.info("Fetching Analysis details.")
        analysis_analysis_alias = self._fetch_analysis_alias(study_accession)
        analysis_analysis_title = "" #TODO: fix this as DGVA's ANALYSIS_ID
        analysis_analysis_description = self._fetch_analysis_description(study_accession)
        analysis_experiment_type = self._fetch_experiment_type(study_accession)
        analysis_reference_genome = self._fetch_reference_genome(study_accession)
        analysis_evidence_type = self._determine_evidence_type(vcf_output)
        analysis_reference_fasta = assembly
        analysis_assembly_report = assembly_report
        analysis_platform = self._fetch_analysis_platform(study_accession)
        analysis_software = self._fetch_analysis_software(study_accession)
        analysis_pipeline_descriptions = self._fetch_analysis_pipeline_descriptions(study_accession)
        analysis_imputation = False # assumption
        analysis_phasing = False # assumption
        analysis_date = "" # not present in DGVA
        analysis_centre = "" # if different to project centre
        analysis_links = self._fetch_analysis_links(study_accession) #TODO: should be [] and  match regex (DB:ID:LABEL).
        analysis_run_accessions = self._fetch_analysis_run_accessions(study_accession)

        # cleaning
        analysis_run_accessions = [r if r else "" for r in analysis_run_accessions]

        # validation
        run_accession_pattern = "^(E|D|S)RR[0-9]{6,}$"
        for run_accession in analysis_run_accessions:
            if run_accession is not None and run_accession != "":
                assert re.fullmatch(run_accession_pattern, run_accession), f"String {run_accession} does not match pattern: {run_accession_pattern}"

        # required: analysisTitle, analysisAlias, description, experimentType, referenceGenome
        # not required: evidenceType, referenceFasta, assemblyReport, platform, software, pipelineDescriptions
        # imputation, phasing, date, links, runAccessions
        analysis_array = []
        analysis_object = {
            "analysisTitle": analysis_analysis_title,
            "analysisAlias": analysis_analysis_alias,
            "description": analysis_analysis_description,
            "experimentType": analysis_experiment_type,
            "referenceGenome": analysis_reference_genome,
            "evidenceType": analysis_evidence_type,
            "referenceFasta": analysis_reference_fasta,
            "assemblyReport": analysis_assembly_report,
            "platform": analysis_platform,
            "software": analysis_software,
            "pipelineDescriptions": analysis_pipeline_descriptions,
            "imputation": analysis_imputation,
            "phasing": analysis_phasing,
            "date": analysis_date,
            "centre": analysis_centre,
            "links":analysis_links,
            "runAccessions": analysis_run_accessions
        }
        placeholder_keys = [k for k,v in analysis_object.items() if v is None or v == ""]
        for placeholder_key in placeholder_keys:
            logger.info(f"{placeholder_key} not found. Adding placeholder: {analysis_object[placeholder_key]}")
        analysis_array.append(analysis_object)
        return analysis_array

    def _get_sample_pre_registered(self, study_accession, biosample_accession, sample_id):
        # requires analysisAlias, sampleinVCF, biosample_accession
        sample_analysis_alias = self._fetch_analysis_alias(study_accession)
        # assumption the name: sampleinVCF = sample_id
        sample_sampleinvcf = sample_id
        sample_object = {
            "analysisAlias": [sample_analysis_alias],
            "sampleInVCF": sample_sampleinvcf,
            "bioSampleAccession": biosample_accession
        }
        return sample_object

    def _get_sample_new(self, study_accession, sample_id):
        # return sample object
        # requires analysisAlias, sampleinVCF, bioSampleObject
        # bioSampleObject requires = name, taxID, scientific_name, release (hold-date) which can be found in DGVA
        # bioSampleObject requires = collection date, geo loc which can be set to unknown/not collected
        sample_analysis_alias_list = self._fetch_analysis_alias_list(study_accession)
        if not sample_analysis_alias_list:
            sample_analysis_alias_list.append("")
        # assume sample in VCF = sample_id
        sample_sampleinvcf = sample_id
        # TODO: fix the error, as there multiple sample names,
        sample_tax_id = self._fetch_tax_id(study_accession)
        scientific_name = self._fetch_scientific_name(study_accession)
        hold_date = self._fetch_hold_date(study_accession)
        collection_date = "not provided"
        geographic_location_country_and_or_sea = "not provided"
        biosample_object = {
            "sample_title": sample_id,
            "scientific_name": scientific_name,
            "tax_id": sample_tax_id,
            "collection date": collection_date,
            "geographic location (country and/or sea)": geographic_location_country_and_or_sea
        }
        sample_object = {
            "analysisAlias": sample_analysis_alias_list,
            "sampleInVCF": sample_sampleinvcf,
            "bioSampleObject": biosample_object
        }
        return sample_object

    def _get_files(self, study_accession, vcf_output):
        #TODO: the files mentioned in DGVA are XML files not the GVF files we are looking for.

        # return files_array
        # requires analysisAlias, fileName

        files_analysis_alias = self._fetch_analysis_alias(study_accession)
        files_file_name = self._get_file_name(vcf_output)
        files_file_size = self._get_file_size(vcf_output)
        files_file_md5 = self._get_file_md5(vcf_output)
        files_array = []
        # required: analysisAlias, filename
        # not required: file size, md5
        files_object = {
            "analysisAlias": files_analysis_alias,
            "fileName": files_file_name,
            "fileSize": files_file_size,
            "md5": files_file_md5
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
        #TODO: change this. Discussion 18/3/26: moving towards treating BIOPROJECT_ACCESSION as a parent project. If treated as pre-registered we will add an analysis to the project.
        #TODO: find how to identify pre-registered projects in DGVA. Will there be pre-registered projects?
        project_accession_query = (Query
                                   .from_(ds)
                                   .select(ds.BIOPROJECT_ACCESSION)
                                   .where(ds.STUDY_ACCESSION == study_accession)
                                   )
        project_accession_dict = self.load_from_db(project_accession_query.get_sql(quote_char=None))
        project_accession = self.validate_fetch_result("projectAccession", project_accession_dict)
        is_project_preregistered = False
        if project_accession is not None:
            is_project_preregistered = True
            logger.info(f"Determining if the project is pre-registered - SUCCESS - Project found: {project_accession}.")
        else:
            logger.info(f"Determining if the project is pre-registered - FAILURE - No project found.")
        return is_project_preregistered, project_accession

    def _determine_sample_pre_registered(self, study_accession):
        # create the schema objects
        db = Schema("DGVA")
        # create the table objects
        dsamp = Table("DGVA_SAMPLE", schema=db).as_("dsamp")
        # create the query
        sample_accession_query = (Query
                                   .from_(dsamp)
                                   .select(dsamp.BIOSAMPLE_ACCESSION, dsamp.SUBMITTER_SAMPLE_ID)
                                   .where(dsamp.STUDY_ACCESSION == study_accession)
                                   )
        # Index: (BiosampleAccession, SubmitterSampleID)
        sample_accession_dict = self.load_from_db(sample_accession_query.get_sql(quote_char=None))
        sample_accession_and_status_list = []
        SampleStatus = namedtuple('SampleStatus', ['is_sample_preregistered', 'sample_accession', 'sample_id'])

        is_sample_preregistered = False
        # key: index, value= Tuple(BiosampleAccession,SampleID) e.g. (None, 'NA20344')
        for value in sample_accession_dict.values():
            if isinstance(value, tuple) and len(value) > 0:
                current_sample_accession = value[0]
                current_sample_id = value[1]
                if current_sample_accession is not None:
                    # DGVa does not have any samples with a BioSample accession. We do not expect this to be populated.
                    is_sample_preregistered = True
                    sample_accession_and_status_list.append(SampleStatus(is_sample_preregistered, current_sample_accession, current_sample_id))
                    # return is_sample_preregistered, sample_accession, sample_id
                    logger.info(f"Determining if sample is pre-registered - SUCCESS - Sample found: {SampleStatus(is_sample_preregistered, current_sample_accession, current_sample_id)}.")
                else:
                    sample_accession_and_status_list.append(SampleStatus(is_sample_preregistered, None, current_sample_id))
                    # return is_sample_preregistered, None, sample_id
                    logger.info(f"Determining if sample is pre-registered - FAILURE - Sample not found: {current_sample_id}.")
        return sample_accession_and_status_list

    def _determine_evidence_type(self, vcf_output):
        with open(vcf_output, "r") as vcf:
            for line in vcf:
                if line.startswith("#CHROM"):
                    header_tokens = line.split("\t")
        number_of_header_tokens = len(header_tokens)
        if number_of_header_tokens == 8:
            evidence_type = "allele_frequency"
        else:
            evidence_type = "genotype"
        logger.info(f"{number_of_header_tokens} tokens found in the VCF header. Determining evidence type as: {evidence_type}")
        return evidence_type

    # Formatting
    def _format_project(self, project_hold_date, project_links, project_parent_project, project_publications):
        if project_hold_date is None:
            project_hold_date = ""
        if project_links != "" and project_links is not None:
            project_links = project_links + "| URL"
        else:
            project_links = ""
        pubmed_publications = []
        if project_parent_project is None:
            project_parent_project = ""
        if type(project_publications) is not list:
            project_publications = [project_publications]
        for pub in project_publications:
            pubmed_string = "PubMed:" + str(pub)
            pubmed_publications.append(pubmed_string)
        return project_hold_date, project_links, project_parent_project, pubmed_publications
    # VALIDATING
    def validate_fetch_result(self, eva_field_name, fetch_result_dict):
        try:
            if fetch_result_dict:
                results = [v[0] for v in fetch_result_dict.values() if v]
                fetch_result = results[0] if len(results) == 1 else results
                if fetch_result:
                    # SUCCESS if value is present or None
                    logger.info(f"Fetching {eva_field_name} - SUCCESS - Value(s) for {eva_field_name} found: {fetch_result}.")
            else:
                raise ValueError(f"Missing data: {eva_field_name}.")
        except ValueError as e:
            logger.error(f"Fetching {eva_field_name}  - FAILURE - {eva_field_name} not found. {e} Setting value as empty string.")
            fetch_result = ""
        return fetch_result
        # try:
        #     if fetch_result_dict:
        #         print(f"The fetch_result_dict for eva_field_name {eva_field_name}: {len(fetch_result_dict)}")
        #         value_list = next(iter(fetch_result_dict.values()), [None])
        #         fetch_result = value_list[0] if value_list else None
        #         if fetch_result:
        #             # SUCCESS if value is present or None
        #             logger.info(f"Fetching {eva_field_name} - SUCCESS - {eva_field_name} found: {fetch_result}.")
        #     else:
        #         raise ValueError(f"Missing data: {eva_field_name}.")
        # except ValueError as e:
        #     logger.error(f"Fetching {eva_field_name}  - FAILURE - {eva_field_name} not found. {e} Setting value as empty string.")
        #     fetch_result = ""
        # return fetch_result

    def validate_date(self, date):
        try:
            datetime.strptime(date, '%Y-%m-%d')
            return True
        except ValueError:
            return False

    def validate_project(self, project_description, project_hold_date, project_parent_project,
                         project_tax_id, project_title, pubmed_publications):
        """ Asserts whether the input parameters meet the EVA JSON schema https://github.com/EBIvariation/eva-sub-cli/blob/main/eva_sub_cli/etc/eva_schema.json
        :params: project_description: string of max 5000 chars
        :params: project_hold_date: YYYY-MM-DD or ""
        :params: project_parent_project: project accession matching regex "^PRJ(E|D|N)[A-Z][0-9]+$"
        :params: project_tax_id: taxonomy ID as an integer
        :params: project_title: string of max 500 chars
        :params: pubmed_publications: list of pubmed ids ['Pubmed:1239234'] (can be more than one, in some studies)
        """
        # constants
        MAX_PROJECT_DESCRIPTION_LENGTH = 5000
        MAX_PROJECT_TITLE_LENGTH = 500
        project_accession_pattern = r"^PRJ(E|D|N)[A-Z][0-9]+$"  # applies to parent/child/peer
        publications_pattern = "^[^:,]+?:[^:,]+?$"  # e.g. PubMed:23128226
        # performing checks
        assert len(
            project_description) <= MAX_PROJECT_DESCRIPTION_LENGTH, f"Project description exceeded length: {MAX_PROJECT_DESCRIPTION_LENGTH}"
        assert len(
            project_title) <= MAX_PROJECT_TITLE_LENGTH, f"Project title exceeded length: {MAX_PROJECT_TITLE_LENGTH}"
        assert isinstance(project_tax_id,
                          int), f"Project Tax ID must be an int: {project_tax_id} is {type(project_tax_id)}"
        if project_hold_date != "":
            assert self.validate_date(
                project_hold_date) == True, f"Project Hold Date must be YYYY-MM-DD: {project_hold_date}"
        if project_parent_project is not None and project_parent_project != "":
            assert re.fullmatch(project_accession_pattern, project_parent_project), f"String {project_parent_project} does not match pattern: {project_accession_pattern}"
        for pub in pubmed_publications:
            assert re.fullmatch(publications_pattern,
                                pub), f"String {pub} does not match pattern: {publications_pattern}"

    #### THESE FETCH FIELDS FROM THE DB
    # SUBMITTER DETAILS SECTION
    def _fetch_submitter_details_dictionary(self, study_accession):
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
            .select(
                sc.FIRST_NAME,
                sc.LAST_NAME,
                sc.CONTACT_PHONE,
                sc.CONTACT_EMAIL,
                sc.AFFILIATION_NAME,
                sc.AFFILIATION_ADDRESS
            ).where(ds.STUDY_ACCESSION == study_accession)
        )
        submitter_details_dict = self.load_from_db(submitter_details_query.get_sql(quote_char=None))
        return submitter_details_dict
    # PROJECT SECTION
    def _fetch_project_title(self, study_accession):
        # create the schema objects
        db = Schema("DGVA")
        # create the table objects
        ds = Table("DGVA_STUDY", schema=db).as_("ds")
        project_title_query = (
            Query.from_(ds)
            .select(ds.DISPLAY_NAME)
            .where(ds.STUDY_ACCESSION == study_accession)

        )
        project_title_dict = self.load_from_db(project_title_query.get_sql(quote_char=None))
        project_title = self.validate_fetch_result("title", project_title_dict)
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
        project_description = self.validate_fetch_result("description", project_description_dict)
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
        tax_id = self.validate_fetch_result("taxId", tax_id_dict)
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
        project_centre = self.validate_fetch_result("centre", project_centre_dict)
        return project_centre

    def _fetch_project_publications(self, study_accession):
        # some study accessions do have multiple publications e.g. estd192
        # create the schema objects
        db = Schema("DGVA")
        # create the table objects
        spp = Table("STUDY_PUBMED_PUBLICATION", schema=db).as_("spp")
        project_publications_query = (
            Query.from_(spp)
            .select(spp.PUBMED_ID)
            .where(spp.STUDY_ACCESSION == study_accession)
        )
        project_publications_dict = self.load_from_db(project_publications_query.get_sql(quote_char=None))
        project_publications = self.validate_fetch_result("publications", project_publications_dict)
        return project_publications

    def _fetch_project_parent_project(self, study_accession):
        # create the schema objects
        db = Schema("DGVA")
        # create the table objects
        ds = Table("DGVA_STUDY", schema=db).as_("ds")
        # create the query
        parent_project_query = (Query
                                   .from_(ds)
                                   .select(ds.BIOPROJECT_ACCESSION)
                                   .where(ds.STUDY_ACCESSION == study_accession)
                                   )
        parent_project_dict = self.load_from_db(parent_project_query.get_sql(quote_char=None))
        parent_project = self.validate_fetch_result("projectAccession", parent_project_dict)
        return parent_project

    def _fetch_project_links(self, study_accession):
        # study url (not experimentlinkurl - that goes in Analysis section)
        # create the schema objects
        db = Schema("DGVA")
        # create the table objects
        ds = Table("DGVA_STUDY", schema=db).as_("ds")
        sub = Table("DGVA_SUBMISSION", schema=db).as_("sub")
        project_links_query = (
            Query.from_(ds)
            .join(sub)
            .on(ds.STUDY_ACCESSION == sub.STUDY_ACCESSION)
            .select(ds.STUDY_URL)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        project_links_dict = self.load_from_db(project_links_query.get_sql(quote_char=None))
        project_links = self.validate_fetch_result("links", project_links_dict)
        return project_links

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
        experiment_type = self.validate_fetch_result("experimentType", experiment_type_dict)
        return experiment_type

    def _fetch_reference_genome(self, study_accession):
        # TODO: reference genome is not in GCA needs converting
        pass

    def _fetch_analysis_platform(self, study_accession):
        # create the schema objects
        db = Schema("DGVA")
        # create the table objects
        de = Table("DGVA_EXPERIMENT", schema=db).as_("de")
        ep = Table("EXPERIMENT_PLATFORM", schema=db).as_("ep")
        dp = Table("DGVA_PLATFORM", schema=db).as_("dp")
        analysis_platform_query = (
            Query.from_(de)
            .join(ep).on(de.EXPERIMENT_ID == ep.EXPERIMENT_ID)
            .join(dp).on(dp.PLATFORM_ID == ep.PLATFORM_ID)
            .select(dp.PLATFORM_NAME)
            .where(de.STUDY_ACCESSION == study_accession)
        )
        analysis_platform_dict = self.load_from_db(analysis_platform_query.get_sql(quote_char=None))
        analysis_platform = self.validate_fetch_result("platform", analysis_platform_dict)
        return analysis_platform

    def _fetch_analysis_software(self, study_accession):
        # create the schema objects
        db = Schema("DGVA")
        # create the table objects
        dd = Table("DGVA_DETECTION", schema=db).as_("dd")
        ed = Table("EXPERIMENT_DETECTION", schema=db).as_("ed")
        de = Table("DGVA_EXPERIMENT", schema=db).as_("de")
        analysis_software_query = (
            Query.from_(dd)
            .join(ed).on(dd.DETECTION_ID == ed.DETECTION_ID)
            .join(de).on(ed.EXPERIMENT_ID == de.EXPERIMENT_ID)
            .select(dd.DETECTION_METHOD)
            .where(de.STUDY_ACCESSION == study_accession)
        )
        analysis_software_dict = self.load_from_db(analysis_software_query.get_sql(quote_char=None))
        analysis_software = self.validate_fetch_result("software", analysis_software_dict)
        return analysis_software

    def _fetch_analysis_pipeline_descriptions(self, study_accession):
        # create the schema objects
        db = Schema("DGVA")
        # create the table objects
        dd = Table("DGVA_DETECTION", schema=db).as_("dd")
        ed = Table("EXPERIMENT_DETECTION", schema=db).as_("ed")
        de = Table("DGVA_EXPERIMENT", schema=db).as_("de")
        analysis_pipeline_descriptions_query = (
            Query.from_(dd)
            .join(ed).on(dd.DETECTION_ID == ed.DETECTION_ID)
            .join(de).on(ed.EXPERIMENT_ID == de.EXPERIMENT_ID)
            .select(dd.DETECTION_DESCRIPTION)
            .where(de.STUDY_ACCESSION == study_accession)
        )
        analysis_pipeline_descriptions_dict = self.load_from_db(analysis_pipeline_descriptions_query.get_sql(quote_char=None))
        analysis_pipeline_descriptions = self.validate_fetch_result("pipelineDescriptions", analysis_pipeline_descriptions_dict)
        return analysis_pipeline_descriptions

    def _fetch_analysis_links(self, study_accession):
        # (SELECT elu.URL FROM DGVA_EXPERIMENT de JOIN EXPERIMENT_LINK_URL elu ON de.EXPERIMENT_ID = elu.EXPERIMENT_ID)
        # create the schema objects
        db = Schema("DGVA")
        # create the table objects
        de = Table("DGVA_EXPERIMENT", schema=db).as_("de")
        elu = Table("EXPERIMENT_LINK_URL", schema=db).as_("elu")
        analysis_links_query = (
            Query.from_(de)
            .join(elu).on(de.EXPERIMENT_ID == elu.EXPERIMENT_ID)
            .select(elu.URL)
            .where(de.STUDY_ACCESSION == study_accession)
        )
        analysis_links_dict = self.load_from_db(analysis_links_query.get_sql(quote_char=None))
        analysis_links = self.validate_fetch_result("links", analysis_links_dict)
        return analysis_links

    def _fetch_analysis_run_accessions(self, study_accession):
        # create the schema objects
        db = Schema("DGVA")
        # create the table objects
        de = Table("DGVA_EXPERIMENT", schema=db).as_("de")
        ds = Table("DGVA_STUDY", schema=db).as_("ds")
        analysis_run_accessions_query = (
            Query.from_(de)
            .join(ds).on(de.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .select(de.SRA_ACCESSION)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        analysis_run_accessions_dict = self.load_from_db(analysis_run_accessions_query.get_sql(quote_char=None))
        analysis_run_accessions = self.validate_fetch_result("runAccessions", analysis_run_accessions_dict)
        return analysis_run_accessions
    # SAMPLE SECTION
    def _fetch_analysis_alias_list(self, study_accession):
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
        analysis_alias_list = [v[0] for v in analysis_alias_dict.values()]
        return analysis_alias_list

    def _fetch_sample_id_list(self, study_accession):
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
        sample_id_list = [v[0] for v in sample_id_dict.values()]
        return sample_id_list

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
        hold_date = self.validate_fetch_result("holdDate", hold_date_dict)
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
        scientific_name = self.validate_fetch_result("scientific_name", scientific_name_dict)
        return scientific_name

    # FILES SECTION
    # THESE GETTERS GET THE RELEVANT VALUE
    def _get_file_name(self, vcf_output):
        file_name = os.path.basename(vcf_output)
        return file_name

    def _get_file_size(self, vcf_output):
        file_size = os.path.getsize(vcf_output)
        return file_size

    def _get_file_md5(self, vcf_output):
        hash_md5 = hashlib.md5()
        with open(vcf_output, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()
