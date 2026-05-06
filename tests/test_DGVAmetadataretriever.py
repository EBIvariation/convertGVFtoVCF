import json
import os
from unittest import TestCase
from unittest.mock import patch, MagicMock, Mock, mock_open

from convert_gvf_to_vcf.metadata_retrievers.dgvametadata import DGVAMetadataRetriever
from convert_gvf_to_vcf.metadata_retrievers.evametadata import EVAMetadataRetriever
from convert_gvf_to_vcf.projectpaths import ProjectPaths

class TestDGVAMetadataRetriever(TestCase):
    def setUp(self):
        input_folder = os.path.dirname(__file__)
        self.config = os.path.join(input_folder, "input", "test.config")
        # self.vcf_output = os.path.join(input_folder, "output", "a.vcf")
        self.json_output = os.path.join(input_folder, "output", "a_dgva.json")
        self.paths = ProjectPaths()
        self.full_config_path = self.paths.full_config_path

    def test_init(self):
        # testing the __init__ function has loaded the config correctly
        metadata_client = DGVAMetadataRetriever(self.config)
        assert metadata_client._host == "mockhost"
        assert metadata_client._port == "1111"
        assert metadata_client._username == "mockuser"
        assert metadata_client._password == "mockpassword"
        assert metadata_client._service_name == "mockservicename"
        assert metadata_client._connection == None
        assert metadata_client._max_retries == 3

    @patch("convert_gvf_to_vcf.metadata_retrievers.dgvametadata.DGVAMetadataRetriever._fetch_method_type")
    @patch("convert_gvf_to_vcf.metadata_retrievers.dgvametadata.DGVAMetadataRetriever._fetch_curated_set_link")
    @patch("convert_gvf_to_vcf.metadata_retrievers.dgvametadata.DGVAMetadataRetriever._fetch_curated_set_name")
    @patch("convert_gvf_to_vcf.metadata_retrievers.dgvametadata.DGVAMetadataRetriever._fetch_curator_email")
    @patch("convert_gvf_to_vcf.metadata_retrievers.dgvametadata.DGVAMetadataRetriever._fetch_curator_name")
    @patch("convert_gvf_to_vcf.metadata_retrievers.dgvametadata.DGVAMetadataRetriever._fetch_detection_description")
    @patch("convert_gvf_to_vcf.metadata_retrievers.dgvametadata.DGVAMetadataRetriever._fetch_detection_method")
    @patch("convert_gvf_to_vcf.metadata_retrievers.dgvametadata.DGVAMetadataRetriever._fetch_experiment_resolution")
    @patch("convert_gvf_to_vcf.metadata_retrievers.dgvametadata.DGVAMetadataRetriever._fetch_experiment_site")
    @patch("convert_gvf_to_vcf.metadata_retrievers.dgvametadata.DGVAMetadataRetriever._fetch_affiliation_url")
    @patch("convert_gvf_to_vcf.metadata_retrievers.dgvametadata.DGVAMetadataRetriever._fetch_correction")
    @patch("convert_gvf_to_vcf.metadata_retrievers.dgvametadata.DGVAMetadataRetriever._fetch_new_feature")
    @patch("convert_gvf_to_vcf.metadata_retrievers.dgvametadata.DGVAMetadataRetriever._fetch_update_comment")
    @patch("convert_gvf_to_vcf.metadata_retrievers.dgvametadata.DGVAMetadataRetriever._fetch_submission_version")
    @patch("convert_gvf_to_vcf.metadata_retrievers.dgvametadata.DGVAMetadataRetriever._fetch_comment_timestamp")
    @patch("convert_gvf_to_vcf.metadata_retrievers.dgvametadata.DGVAMetadataRetriever._fetch_comment_user_name")
    @patch("convert_gvf_to_vcf.metadata_retrievers.dgvametadata.DGVAMetadataRetriever._fetch_study_comment")
    @patch("convert_gvf_to_vcf.metadata_retrievers.dgvametadata.DGVAMetadataRetriever._fetch_creation_date")
    def test_create_json_dgva(self,
                              mock_creation_date,
                              mock_study_comment,
                              mock_comment_user_name,
                              mock_comment_timestamp,
                              mock_submission_version,
                              mock_update_comment,
                              mock_new_feature,
                              mock_correction,
                              mock_affiliation_url,
                              mock_experiment_site,
                              mock_experiment_resolution,
                              mock_detection_method,
                              mock_detection_description,
                              mock_curator_name,
                              mock_curation_email,
                              mock_curated_set_name,
                              mock_curated_set_link,
                              mock_method_type
                              ):
        mock_map = {
                              mock_creation_date: "mock_creation_date",
                              mock_study_comment: "mock_study_comment",
                              mock_comment_user_name: "mock_comment_user_name",
                              mock_comment_timestamp: "mock_comment_timestamp",
                              mock_submission_version: 1,
                              mock_update_comment: "mock_update_comment",
                              mock_new_feature: "mock_new_feature",
                              mock_correction: 0,
                              mock_affiliation_url: "mock_affiliation_url",
                              mock_experiment_site: "mock_experiment_site",
                              mock_experiment_resolution: "mock_experiment_resolution",
                              mock_detection_method: "mock_detection_method",
                              mock_detection_description: "mock_detection_description",
                              mock_curator_name: "mock_curator_name",
                              mock_curation_email: "mock_curation_email",
                              mock_curated_set_name: "mock_curated_set_name",
                              mock_curated_set_link: "mock_curated_set_link",
                              mock_method_type: "mock_method_type"
        }

        # set up mock result
        for mock, label in mock_map.items():
            mock.return_value = f"testing_{label}"

        metadata_client = DGVAMetadataRetriever(self.config)
        metadata_client.create_json_dgva(self.json_output, study_accession='estd22')

        for mock in mock_map.keys():
            mock.assert_called_once()
        assert os.path.exists(self.json_output)
        with open(self.json_output, 'r') as f:
            results = json.load(f)

        expected = {'dgva': [{'creationDate': 'testing_mock_creation_date', 'studyComment': 'testing_mock_study_comment', 'commentUserName': 'testing_mock_comment_user_name', 'commentTimestamp': 'testing_mock_comment_timestamp', 'submissionVersion': 'testing_mock_submission_version', 'updateComment': 'testing_mock_update_comment', 'newFeature': 'testing_mock_new_feature', 'correction': 'testing_mock_correction', 'affiliationUrl': 'testing_mock_affiliation_url', 'experimentSite': 'testing_mock_experiment_site', 'experimentResolution': 'testing_mock_experiment_resolution', 'detectionMethod': 'testing_mock_detection_method', 'detectionDescription': 'testing_mock_detection_description', 'curatorName': 'testing_mock_curator_name', 'curatorEmail': 'testing_mock_curation_email', 'curatedSetName': 'testing_mock_curated_set_name', 'curatedSetLink': 'testing_mock_curated_set_link', 'methodType': 'testing_mock_method_type'}]}
        self.assertEqual(results, expected)

