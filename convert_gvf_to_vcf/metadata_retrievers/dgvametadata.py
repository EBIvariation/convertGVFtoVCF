import json
from jsonschema import validate, ValidationError
from pypika import Query, Table, Schema

from ebi_eva_common_pyutils.logger import logging_config as log_cfg

logger = log_cfg.get_logger(__name__)

from convert_gvf_to_vcf.metadata_retrievers.basemetadataretriever import BaseMetadataRetriever
class DGVAMetadataRetriever(BaseMetadataRetriever):
    def retrieve(self):
        pass

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
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        project_title_query = (
            Query.from_(ds)
            .select(ds.CREATION_DATE)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        creation_date_list = self.load_from_db(project_title_query.get_sql(quote_char=None))
        [creation_date, *_] = self.fetch_results_from_rows("creationDate", creation_date_list) or [""]
        creation_date_json_string = json.dumps(creation_date, default=str)
        return creation_date_json_string or ""

    def _fetch_study_comment(self, study_accession):
        # create the table objects
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        sc = Table("STUDY_COMMENT", schema=self.db).as_("sc")
        study_comment_query=(
            Query.from_(sc)
            .join(ds)
            .on(sc.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .select(sc.STUDY_COMMENT)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        study_comment_list = self.load_from_db(study_comment_query.get_sql(quote_char=None))
        [study_comment, *_] = self.fetch_results_from_rows("studyComment", study_comment_list) or [""]
        return study_comment or ""

    def _fetch_comment_user_name(self, study_accession):
        # create the table objects
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        sc = Table("STUDY_COMMENT", schema=self.db).as_("sc")
        study_comment_query=(
            Query.from_(sc)
            .join(ds)
            .on(sc.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .select(sc.COMMENT_USER_NAME)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        comment_user_name_list = self.load_from_db(study_comment_query.get_sql(quote_char=None))
        [comment_user_name, *_] = self.fetch_results_from_rows("commentUserName", comment_user_name_list) or [""]
        return comment_user_name or ""

    def _fetch_comment_timestamp(self, study_accession):
        # create the table objects
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        sc = Table("STUDY_COMMENT", schema=self.db).as_("sc")
        study_comment_query=(
            Query.from_(sc)
            .join(ds)
            .on(sc.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .select(sc.COMMENT_TIMESTAMP)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        comment_timestamp_list = self.load_from_db(study_comment_query.get_sql(quote_char=None))
        [comment_timestamp, *_] = self.fetch_results_from_rows("commentTimestamp", comment_timestamp_list) or [""]
        return comment_timestamp or ""

    def _fetch_submission_version(self, study_accession):
        # create the table objects
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        su = Table("STUDY_UPDATE", schema=self.db).as_("su")
        study_comment_query=(
            Query.from_(su)
            .join(ds)
            .on(su.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .select(su.SUBMISSION_VERSION)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        submission_version_list = self.load_from_db(study_comment_query.get_sql(quote_char=None))
        [submission_version, *_] = self.fetch_results_from_rows("submissionVersion", submission_version_list) or [""]
        return submission_version or ""

    def _fetch_update_comment(self, study_accession):
        # create the table objects
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        su = Table("STUDY_UPDATE", schema=self.db).as_("su")
        study_comment_query=(
            Query.from_(su)
            .join(ds)
            .on(su.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .select(su.UPDATE_COMMENT)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        update_comment_list = self.load_from_db(study_comment_query.get_sql(quote_char=None))
        [update_comment, *_] = self.fetch_results_from_rows("updateComment", update_comment_list) or [""]
        return update_comment or ""

    def _fetch_new_feature(self, study_accession):
        # create the table objects
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        su = Table("STUDY_UPDATE", schema=self.db).as_("su")
        study_comment_query=(
            Query.from_(su)
            .join(ds)
            .on(su.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .select(su.NEW_FEATURE)
            .where(ds.STUDY_ACCESSION == study_accession)
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
        de = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        dm = Table("DGVA_METHOD", schema=self.db).as_("dm")
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        em = Table("EXPERIMENT_METHOD", schema=self.db).as_("em")
        method_type_query=(
            Query.from_(dm)
            .join(em)
            .on(dm.METHOD_ID == em.METHOD_ID)
            .join(de)
            .on(de.EXPERIMENT_ID == em.EXPERIMENT_ID)
            .join(ds)
            .on(dm.METHOD_ID == em.METHOD_ID)
            .select(dm.METHOD_TYPE)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        method_type_list = self.load_from_db(method_type_query.get_sql(quote_char=None))
        [method_type, *_] = self.fetch_results_from_rows("methodType", method_type_list) or [""]
        return method_type or ""

    def _fetch_curated_set_link(self, study_accession):
        # create the table objects
        de = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        ec = Table("EXPERIMENT_CURATION", schema=self.db).as_("ec")
        curated_set_link_query=(
            Query.from_(ec)
            .join(de)
            .on(de.EXPERIMENT_ID == ec.EXPERIMENT_ID)
            .join(ds)
            .on(ds.STUDY_ACCESSION == de.STUDY_ACCESSION)
            .select(ec.CURATED_SET_LINK)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        curated_set_link_list = self.load_from_db(curated_set_link_query.get_sql(quote_char=None))
        [curated_set_link, *_] = self.fetch_results_from_rows("curatedSetLink", curated_set_link_list) or [""]
        return curated_set_link or ""

    def _fetch_curated_set_name(self, study_accession):
        # create the table objects
        de = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        ec = Table("EXPERIMENT_CURATION", schema=self.db).as_("ec")
        curated_set_name_query=(
            Query.from_(ec)
            .join(de)
            .on(de.EXPERIMENT_ID == ec.EXPERIMENT_ID)
            .join(ds)
            .on(ds.STUDY_ACCESSION == de.STUDY_ACCESSION)
            .select(ec.CURATED_SET_NAME)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        curated_set_name_list = self.load_from_db(curated_set_name_query.get_sql(quote_char=None))
        [curated_set_name, *_] = self.fetch_results_from_rows("curatedSetName", curated_set_name_list) or [""]
        return curated_set_name or ""

    def _fetch_curator_email(self, study_accession):
        # create the table objects
        de = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        ec = Table("EXPERIMENT_CURATION", schema=self.db).as_("ec")
        curator_email_query = (
            Query.from_(ec)
            .join(de)
            .on(de.EXPERIMENT_ID == ec.EXPERIMENT_ID)
            .join(ds)
            .on(ds.STUDY_ACCESSION == de.STUDY_ACCESSION)
            .select(ec.CURATOR_EMAIL)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        curator_email_list = self.load_from_db(curator_email_query.get_sql(quote_char=None))
        [curator_email, *_] = self.fetch_results_from_rows("curatorEmailList", curator_email_list) or [""]
        return curator_email or ""

    def _fetch_curator_name(self, study_accession):
        # create the table objects
        de = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        ec = Table("EXPERIMENT_CURATION", schema=self.db).as_("ec")
        curator_name_query = (
            Query.from_(ec)
            .join(de)
            .on(de.EXPERIMENT_ID == ec.EXPERIMENT_ID)
            .join(ds)
            .on(ds.STUDY_ACCESSION == de.STUDY_ACCESSION)
            .select(ec.CURATOR_NAME)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        curator_name_list = self.load_from_db(curator_name_query.get_sql(quote_char=None))
        [curator_name, *_] = self.fetch_results_from_rows("curatorNameList", curator_name_list) or [""]
        return curator_name or ""

    def _fetch_detection_description(self, study_accession):
        # create the table objects
        de = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        dd = Table("DGVA_DETECTION", schema=self.db).as_("dd")
        ed = Table("EXPERIMENT_DETECTION", schema=self.db).as_("ed")
        detection_description_query = (
            Query.from_(dd)
            .join(ed)
            .on(dd.DETECTION_ID == ed.DETECTION_ID)
            .join(de)
            .on(de.EXPERIMENT_ID == ed.EXPERIMENT_ID)
            .join(ds)
            .on(ds.STUDY_ACCESSION == de.STUDY_ACCESSION)
            .select(dd.DETECTION_DESCRIPTION)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        detection_description_list = self.load_from_db(detection_description_query.get_sql(quote_char=None))
        [detection_description, *_] = self.fetch_results_from_rows("detectionDescription", detection_description_list) or [""]
        return detection_description or ""

    def _fetch_detection_method(self, study_accession):
        # create the table objects
        de = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        dd = Table("DGVA_DETECTION", schema=self.db).as_("dd")
        ed = Table("EXPERIMENT_DETECTION", schema=self.db).as_("ed")
        detection_method_query = (
            Query.from_(dd)
            .join(ed)
            .on(dd.DETECTION_ID == ed.DETECTION_ID)
            .join(de)
            .on(de.EXPERIMENT_ID == ed.EXPERIMENT_ID)
            .join(ds)
            .on(ds.STUDY_ACCESSION == de.STUDY_ACCESSION)
            .select(dd.DETECTION_METHOD)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        detection_method_list = self.load_from_db(detection_method_query.get_sql(quote_char=None))
        [detection_method, *_] = self.fetch_results_from_rows("detection_method", detection_method_list) or [""]
        return detection_method or ""

    def _fetch_experiment_resolution(self, study_accession):
        # create the table objects
        de = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        experiment_resolution_query = (
            Query.from_(de)
            .join(ds)
            .on(ds.STUDY_ACCESSION == de.STUDY_ACCESSION)
            .select(de.EXPERIMENT_RESOLUTION)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        experiment_resolution_list = self.load_from_db(experiment_resolution_query.get_sql(quote_char=None))
        [experiment_resolution, *_] = self.fetch_results_from_rows("experimentResolution", experiment_resolution_list) or [""]
        return experiment_resolution or ""

    def _fetch_experiment_site(self, study_accession):
        # create the table objects
        de = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        experiment_site_query = (
            Query.from_(de)
            .join(ds)
            .on(ds.STUDY_ACCESSION == de.STUDY_ACCESSION)
            .select(de.EXPERIMENT_SITE)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        experiment_site_list = self.load_from_db(experiment_site_query.get_sql(quote_char=None))
        [experiment_site, *_] = self.fetch_results_from_rows("experimentSite", experiment_site_list) or [""]
        return experiment_site or ""

    def _fetch_affiliation_url(self, study_accession):
        # create the table objects
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        sc = Table("STUDY_CONTACT", schema=self.db).as_("sc")
        study_comment_query=(
            Query.from_(sc)
            .join(ds)
            .on(sc.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .select(sc.AFFILIATION_URL)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        affiliation_url_list = self.load_from_db(study_comment_query.get_sql(quote_char=None))
        [affiliation_url, *_] = self.fetch_results_from_rows("affiliation_url", affiliation_url_list) or [""]
        return affiliation_url or ""

    def _fetch_correction(self, study_accession):
        # create the table objects
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        su = Table("STUDY_UPDATE", schema=self.db).as_("su")
        study_comment_query=(
            Query.from_(su)
            .join(ds)
            .on(su.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .select(su.CORRECTION)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        correction_list = self.load_from_db(study_comment_query.get_sql(quote_char=None))
        [correction, *_] = self.fetch_results_from_rows("correction", correction_list) or [""]
        return correction or 0
