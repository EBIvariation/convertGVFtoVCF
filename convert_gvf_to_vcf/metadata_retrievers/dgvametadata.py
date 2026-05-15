import json
from jsonschema import validate, ValidationError
from pypika import Query, Table, Schema

from ebi_eva_common_pyutils.logger import logging_config as log_cfg

logger = log_cfg.get_logger(__name__)

from convert_gvf_to_vcf.metadata_retrievers.basemetadataretriever import BaseMetadataRetriever
class DGVAMetadataRetriever(BaseMetadataRetriever):

    def create_json_dgva(self, json_file_path, study_accession):
        json_in_dgva_format = {
            "dgva": [
                {
                    "creationDate": self._fetch_creation_date(study_accession),
                    "studyComment": self._fetch_study_comment(study_accession),
                    "commentUserName": self._fetch_comment_user_name(study_accession),
                    "commentTimestamp": self._fetch_comment_timestamp(study_accession),
                    "submissionVersion": self._fetch_submission_version(study_accession),
                    "updateComment": self._fetch_update_comment(study_accession),
                    "newFeature": self._fetch_new_feature(study_accession),
                    "correction": self._fetch_correction(study_accession),
                    "affiliationUrl": self._fetch_affiliation_url(study_accession),
                    "experimentSite": self._fetch_experiment_site(study_accession),
                    "experimentResolution": self._fetch_experiment_resolution(study_accession),
                    "detectionMethod": self._fetch_detection_method(study_accession),
                    "detectionDescription": self._fetch_detection_description(study_accession),
                    "curatorName": self._fetch_curator_name(study_accession),
                    "curatorEmail": self._fetch_curator_email(study_accession),
                    "curatedSetName": self._fetch_curated_set_name(study_accession),
                    "curatedSetLink": self._fetch_curated_set_link(study_accession),
                    "methodType": self._fetch_method_type(study_accession)
                }
            ]
        }
        with open(self.paths.dgva_schema_path, 'r') as dgva_schema_file:
            dgva_schema = json.load(dgva_schema_file)
        validate(instance=json_in_dgva_format, schema=dgva_schema)
        logger.info(f"Validating DGVA JSON file for {study_accession} - SUCCESS")
        with open(json_file_path, 'w') as f:
            json.dump(json_in_dgva_format, f, indent=4)
        logger.info(f"Write DGVA JSON file for {study_accession}- SUCCESS: {json_file_path}")

    def _fetch_creation_date(self, study_accession):
        # create the table objects
        dgva_study_table = Table("DGVA_STUDY", schema=self.db).as_("ds")
        project_title_query = (
            Query.from_(dgva_study_table)
            .select(dgva_study_table.CREATION_DATE)
            .where(dgva_study_table.STUDY_ACCESSION == study_accession)
        )
        creation_date_list = self.load_from_db(project_title_query.get_sql(quote_char=None))
        [creation_date, *_] = self.fetch_results_from_rows("creationDate", creation_date_list) or [""]
        return str(creation_date) if creation_date else ""

    def _fetch_study_comment(self, study_accession):
        # create the table objects
        dgva_study_table = Table("DGVA_STUDY", schema=self.db).as_("ds")
        study_comment_table = Table("STUDY_COMMENT", schema=self.db).as_("sc")
        study_comment_query=(
            Query.from_(study_comment_table)
            .join(dgva_study_table)
            .on(study_comment_table.STUDY_ACCESSION == dgva_study_table.STUDY_ACCESSION)
            .select(study_comment_table.STUDY_COMMENT)
            .where(dgva_study_table.STUDY_ACCESSION == study_accession)
        )
        study_comment_list = self.load_from_db(study_comment_query.get_sql(quote_char=None))
        [study_comment, *_] = self.fetch_results_from_rows("studyComment", study_comment_list) or [""]
        return study_comment or ""

    def _fetch_comment_user_name(self, study_accession):
        # create the table objects
        dgva_study_table = Table("DGVA_STUDY", schema=self.db).as_("ds")
        study_comment_table = Table("STUDY_COMMENT", schema=self.db).as_("sc")
        study_comment_query=(
            Query.from_(study_comment_table)
            .join(dgva_study_table)
            .on(study_comment_table.STUDY_ACCESSION == dgva_study_table.STUDY_ACCESSION)
            .select(study_comment_table.COMMENT_USER_NAME)
            .where(dgva_study_table.STUDY_ACCESSION == study_accession)
        )
        comment_user_name_list = self.load_from_db(study_comment_query.get_sql(quote_char=None))
        [comment_user_name, *_] = self.fetch_results_from_rows("commentUserName", comment_user_name_list) or [""]
        return comment_user_name or ""

    def _fetch_comment_timestamp(self, study_accession):
        # create the table objects
        dgva_study_table = Table("DGVA_STUDY", schema=self.db).as_("ds")
        study_comment_table = Table("STUDY_COMMENT", schema=self.db).as_("sc")
        study_comment_query=(
            Query.from_(study_comment_table)
            .join(dgva_study_table)
            .on(study_comment_table.STUDY_ACCESSION == dgva_study_table.STUDY_ACCESSION)
            .select(study_comment_table.COMMENT_TIMESTAMP)
            .where(dgva_study_table.STUDY_ACCESSION == study_accession)
        )
        comment_timestamp_list = self.load_from_db(study_comment_query.get_sql(quote_char=None))
        [comment_timestamp, *_] = self.fetch_results_from_rows("commentTimestamp", comment_timestamp_list) or [""]
        return str(comment_timestamp) if comment_timestamp else ""

    def _fetch_submission_version(self, study_accession):
        # create the table objects
        dgva_study_table = Table("DGVA_STUDY", schema=self.db).as_("ds")
        study_update_table = Table("STUDY_UPDATE", schema=self.db).as_("su")
        study_comment_query=(
            Query.from_(study_update_table)
            .join(dgva_study_table)
            .on(study_update_table.STUDY_ACCESSION == dgva_study_table.STUDY_ACCESSION)
            .select(study_update_table.SUBMISSION_VERSION)
            .where(dgva_study_table.STUDY_ACCESSION == study_accession)
        )
        submission_version_list = self.load_from_db(study_comment_query.get_sql(quote_char=None))
        [submission_version, *_] = self.fetch_results_from_rows("submissionVersion", submission_version_list) or [""]
        return submission_version or 0

    def _fetch_update_comment(self, study_accession):
        # create the table objects
        dgva_study_table = Table("DGVA_STUDY", schema=self.db).as_("ds")
        study_update_table = Table("STUDY_UPDATE", schema=self.db).as_("su")
        study_comment_query=(
            Query.from_(study_update_table)
            .join(dgva_study_table)
            .on(study_update_table.STUDY_ACCESSION == dgva_study_table.STUDY_ACCESSION)
            .select(study_update_table.UPDATE_COMMENT)
            .where(dgva_study_table.STUDY_ACCESSION == study_accession)
        )
        update_comment_list = self.load_from_db(study_comment_query.get_sql(quote_char=None))
        [update_comment, *_] = self.fetch_results_from_rows("updateComment", update_comment_list) or [""]
        return update_comment or ""

    def _fetch_new_feature(self, study_accession):
        # create the table objects
        dgva_study_table = Table("DGVA_STUDY", schema=self.db).as_("ds")
        study_update_table = Table("STUDY_UPDATE", schema=self.db).as_("su")
        study_comment_query=(
            Query.from_(study_update_table)
            .join(dgva_study_table)
            .on(study_update_table.STUDY_ACCESSION == dgva_study_table.STUDY_ACCESSION)
            .select(study_update_table.NEW_FEATURE)
            .where(dgva_study_table.STUDY_ACCESSION == study_accession)
        )
        new_feature_list = self.load_from_db(study_comment_query.get_sql(quote_char=None))
        [new_feature, *_] = self.fetch_results_from_rows("newFeature", new_feature_list) or [""]
        if new_feature == "1":
            new_feature = True
        else:
            new_feature = False
        return new_feature

    def _fetch_method_type(self, study_accession):
        # create the table objects
        dgva_experiment_table = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        dgva_method_table = Table("DGVA_METHOD", schema=self.db).as_("dm")
        dgva_study_table = Table("DGVA_STUDY", schema=self.db).as_("ds")
        experiment_method_table = Table("EXPERIMENT_METHOD", schema=self.db).as_("em")
        method_type_query=(
            Query.from_(dgva_method_table)
            .join(experiment_method_table)
            .on(dgva_method_table.METHOD_ID == experiment_method_table.METHOD_ID)
            .join(dgva_experiment_table)
            .on(dgva_experiment_table.EXPERIMENT_ID == experiment_method_table.EXPERIMENT_ID)
            .join(dgva_study_table)
            .on(dgva_study_table.STUDY_ACCESSION == dgva_experiment_table.STUDY_ACCESSION)
            .select(dgva_method_table.METHOD_TYPE)
            .where(dgva_study_table.STUDY_ACCESSION == study_accession)
        )
        method_type_list = self.load_from_db(method_type_query.get_sql(quote_char=None))
        [method_type, *_] = self.fetch_results_from_rows("methodType", method_type_list) or [""]
        return method_type or ""

    def _fetch_curated_set_link(self, study_accession):
        # create the table objects
        dgva_experiment_table = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        dgva_study_table = Table("DGVA_STUDY", schema=self.db).as_("ds")
        experiment_curation_table = Table("EXPERIMENT_CURATION", schema=self.db).as_("ec")
        curated_set_link_query=(
            Query.from_(experiment_curation_table)
            .join(dgva_experiment_table)
            .on(dgva_experiment_table.EXPERIMENT_ID == experiment_curation_table.EXPERIMENT_ID)
            .join(dgva_study_table)
            .on(dgva_study_table.STUDY_ACCESSION == dgva_experiment_table.STUDY_ACCESSION)
            .select(experiment_curation_table.CURATED_SET_LINK)
            .where(dgva_study_table.STUDY_ACCESSION == study_accession)
        )
        curated_set_link_list = self.load_from_db(curated_set_link_query.get_sql(quote_char=None))
        [curated_set_link, *_] = self.fetch_results_from_rows("curatedSetLink", curated_set_link_list) or [""]
        return curated_set_link or ""

    def _fetch_curated_set_name(self, study_accession):
        # create the table objects
        dgva_experiment_table = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        dgva_study_table = Table("DGVA_STUDY", schema=self.db).as_("ds")
        experiment_curation_table = Table("EXPERIMENT_CURATION", schema=self.db).as_("ec")
        curated_set_name_query=(
            Query.from_(experiment_curation_table)
            .join(dgva_experiment_table)
            .on(dgva_experiment_table.EXPERIMENT_ID == experiment_curation_table.EXPERIMENT_ID)
            .join(dgva_study_table)
            .on(dgva_study_table.STUDY_ACCESSION == dgva_experiment_table.STUDY_ACCESSION)
            .select(experiment_curation_table.CURATED_SET_NAME)
            .where(dgva_study_table.STUDY_ACCESSION == study_accession)
        )
        curated_set_name_list = self.load_from_db(curated_set_name_query.get_sql(quote_char=None))
        [curated_set_name, *_] = self.fetch_results_from_rows("curatedSetName", curated_set_name_list) or [""]
        return curated_set_name or ""

    def _fetch_curator_email(self, study_accession):
        # create the table objects
        dgva_experiment_table = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        dgva_study_table = Table("DGVA_STUDY", schema=self.db).as_("ds")
        experiment_curation_table = Table("EXPERIMENT_CURATION", schema=self.db).as_("ec")
        curator_email_query = (
            Query.from_(experiment_curation_table)
            .join(dgva_experiment_table)
            .on(dgva_experiment_table.EXPERIMENT_ID == experiment_curation_table.EXPERIMENT_ID)
            .join(dgva_study_table)
            .on(dgva_study_table.STUDY_ACCESSION == dgva_experiment_table.STUDY_ACCESSION)
            .select(experiment_curation_table.CURATOR_EMAIL)
            .where(dgva_study_table.STUDY_ACCESSION == study_accession)
        )
        curator_email_list = self.load_from_db(curator_email_query.get_sql(quote_char=None))
        [curator_email, *_] = self.fetch_results_from_rows("curatorEmailList", curator_email_list) or [""]
        return curator_email or ""

    def _fetch_curator_name(self, study_accession):
        # create the table objects
        dgva_experiment_table = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        dgva_study_table = Table("DGVA_STUDY", schema=self.db).as_("ds")
        experiment_curation_table = Table("EXPERIMENT_CURATION", schema=self.db).as_("ec")
        curator_name_query = (
            Query.from_(experiment_curation_table)
            .join(dgva_experiment_table)
            .on(dgva_experiment_table.EXPERIMENT_ID == experiment_curation_table.EXPERIMENT_ID)
            .join(dgva_study_table)
            .on(dgva_study_table.STUDY_ACCESSION == dgva_experiment_table.STUDY_ACCESSION)
            .select(experiment_curation_table.CURATOR_NAME)
            .where(dgva_study_table.STUDY_ACCESSION == study_accession)
        )
        curator_name_list = self.load_from_db(curator_name_query.get_sql(quote_char=None))
        [curator_name, *_] = self.fetch_results_from_rows("curatorNameList", curator_name_list) or [""]
        return curator_name or ""

    def _fetch_detection_description(self, study_accession):
        # create the table objects
        dgva_experiment_table = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        dgva_study_table = Table("DGVA_STUDY", schema=self.db).as_("ds")
        dgva_detection_table = Table("DGVA_DETECTION", schema=self.db).as_("dd")
        experiment_detection_table = Table("EXPERIMENT_DETECTION", schema=self.db).as_("ed")
        detection_description_query = (
            Query.from_(dgva_detection_table)
            .join(experiment_detection_table)
            .on(dgva_detection_table.DETECTION_ID == experiment_detection_table.DETECTION_ID)
            .join(dgva_experiment_table)
            .on(dgva_experiment_table.EXPERIMENT_ID == experiment_detection_table.EXPERIMENT_ID)
            .join(dgva_study_table)
            .on(dgva_study_table.STUDY_ACCESSION == dgva_experiment_table.STUDY_ACCESSION)
            .select(dgva_detection_table.DETECTION_DESCRIPTION)
            .where(dgva_study_table.STUDY_ACCESSION == study_accession)
        )
        detection_description_list = self.load_from_db(detection_description_query.get_sql(quote_char=None))
        [detection_description, *_] = self.fetch_results_from_rows("detectionDescription", detection_description_list) or [""]
        return detection_description or ""

    def _fetch_detection_method(self, study_accession):
        # create the table objects
        dgva_experiment_table = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        dgva_study_table = Table("DGVA_STUDY", schema=self.db).as_("ds")
        dgva_detection_table = Table("DGVA_DETECTION", schema=self.db).as_("dd")
        experiment_detection_table = Table("EXPERIMENT_DETECTION", schema=self.db).as_("ed")
        detection_method_query = (
            Query.from_(dgva_detection_table)
            .join(experiment_detection_table)
            .on(dgva_detection_table.DETECTION_ID == experiment_detection_table.DETECTION_ID)
            .join(dgva_experiment_table)
            .on(dgva_experiment_table.EXPERIMENT_ID == experiment_detection_table.EXPERIMENT_ID)
            .join(dgva_study_table)
            .on(dgva_study_table.STUDY_ACCESSION == dgva_experiment_table.STUDY_ACCESSION)
            .select(dgva_detection_table.DETECTION_METHOD)
            .where(dgva_study_table.STUDY_ACCESSION == study_accession)
        )
        detection_method_list = self.load_from_db(detection_method_query.get_sql(quote_char=None))
        [detection_method, *_] = self.fetch_results_from_rows("detection_method", detection_method_list) or [""]
        return detection_method or ""

    def _fetch_experiment_resolution(self, study_accession):
        # create the table objects
        dgva_experiment_table = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        dgva_study_table = Table("DGVA_STUDY", schema=self.db).as_("ds")
        experiment_resolution_query = (
            Query.from_(dgva_experiment_table)
            .join(dgva_study_table)
            .on(dgva_study_table.STUDY_ACCESSION == dgva_experiment_table.STUDY_ACCESSION)
            .select(dgva_experiment_table.EXPERIMENT_RESOLUTION)
            .where(dgva_study_table.STUDY_ACCESSION == study_accession)
        )
        experiment_resolution_list = self.load_from_db(experiment_resolution_query.get_sql(quote_char=None))
        [experiment_resolution, *_] = self.fetch_results_from_rows("experimentResolution", experiment_resolution_list) or [""]
        return experiment_resolution or ""

    def _fetch_experiment_site(self, study_accession):
        # create the table objects
        dgva_experiment_table = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        dgva_study_table = Table("DGVA_STUDY", schema=self.db).as_("ds")
        experiment_site_query = (
            Query.from_(dgva_experiment_table)
            .join(dgva_study_table)
            .on(dgva_study_table.STUDY_ACCESSION == dgva_experiment_table.STUDY_ACCESSION)
            .select(dgva_experiment_table.EXPERIMENT_SITE)
            .where(dgva_study_table.STUDY_ACCESSION == study_accession)
        )
        experiment_site_list = self.load_from_db(experiment_site_query.get_sql(quote_char=None))
        [experiment_site, *_] = self.fetch_results_from_rows("experimentSite", experiment_site_list) or [""]
        return experiment_site or ""

    def _fetch_affiliation_url(self, study_accession):
        # create the table objects
        dgva_study_table = Table("DGVA_STUDY", schema=self.db).as_("ds")
        study_contact_table = Table("STUDY_CONTACT", schema=self.db).as_("sc")
        study_comment_query=(
            Query.from_(study_contact_table)
            .join(dgva_study_table)
            .on(study_contact_table.STUDY_ACCESSION == dgva_study_table.STUDY_ACCESSION)
            .select(study_contact_table.AFFILIATION_URL)
            .where(dgva_study_table.STUDY_ACCESSION == study_accession)
        )
        affiliation_url_list = self.load_from_db(study_comment_query.get_sql(quote_char=None))
        [affiliation_url, *_] = self.fetch_results_from_rows("affiliation_url", affiliation_url_list) or [""]
        return affiliation_url or ""

    def _fetch_correction(self, study_accession):
        # create the table objects
        dgva_study_table = Table("DGVA_STUDY", schema=self.db).as_("ds")
        study_update_table = Table("STUDY_UPDATE", schema=self.db).as_("su")
        study_comment_query=(
            Query.from_(study_update_table)
            .join(dgva_study_table)
            .on(study_update_table.STUDY_ACCESSION == dgva_study_table.STUDY_ACCESSION)
            .select(study_update_table.CORRECTION)
            .where(dgva_study_table.STUDY_ACCESSION == study_accession)
        )
        correction_list = self.load_from_db(study_comment_query.get_sql(quote_char=None))
        [correction, *_] = self.fetch_results_from_rows("correction", correction_list) or [""]
        return correction or 0
