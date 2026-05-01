import json
import os
import re
from datetime import datetime

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
                    "creationDate": self._fetch_creation_date(study_accession=study_accession),
                    "studyComment": self._get_study_comment(),
                    "commentUserName": self._get_comment_user_name(),
                    "commentTimestamp": self._get_comment_timestamp(),
                    "submissionVersion": self._get_submission_version(),
                    "updateComment": self._get_update_comment(),
                    "newFeature": self._get_new_feature(),
                    "correction": self._get_correction(),
                    "affiliationUrl": self._get_affiliation_url(),
                    "experimentSite": self._get_experiment_site(),
                    "experimentResolution": self._get_experiment_resolution(),
                    "detectionMethod": self._get_detection_method(),
                    "detectionDescription": self._get_detection_description(),
                    "curatorName": self._get_curator_name(),
                    "curatorEmail": self._get_curator_email(),
                    "curatedSetName": self._get_curated_set_name(),
                    "curatedSetLink": self._get_curated_set_link(),
                    "methodType": self._get_method_type()
                }
            ]
        }
        with open(json_file_path, 'w') as f:
            json.dump(json_in_dgva_format, f, indent=4)
        logger.info(f"Write DGVA JSON file for {study_accession}- SUCCESS: {json_file_path}")

    def _fetch_creation_date(self, study_accession):
        # create the schema objects
        db = Schema("DGVA")
        # create the table objects
        ds = Table("DGVA_STUDY", schema=db).as_("ds")
        project_title_query = (
            Query.from_(ds)
            .select(ds.CREATION_DATE)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        creation_date_list = self.load_from_db(project_title_query.get_sql(quote_char=None))
        [creation_date, *_] = self.fetch_results_from_rows("creationDate", creation_date_list) or [""]
        #TODO: return datetime object in a way json can print
        creation_date = "TESTING"
        return creation_date

    def _get_study_comment(self):
        pass

    def _get_comment_user_name(self):
        pass

    def _get_comment_timestamp(self):
        pass

    def _get_submission_version(self):
        pass

    def _get_update_comment(self):
        pass

    def _get_new_feature(self):
        pass

    def _get_method_type(self):
        pass

    def _get_curated_set_link(self):
        pass

    def _get_curated_set_name(self):
        pass

    def _get_curator_email(self):
        pass

    def _get_curator_name(self):
        pass

    def _get_detection_description(self):
        pass

    def _get_detection_method(self):
        pass

    def _get_experiment_resolution(self):
        pass

    def _get_experiment_site(self):
        pass

    def _get_affiliation_url(self):
        pass

    def _get_correction(self):
        pass