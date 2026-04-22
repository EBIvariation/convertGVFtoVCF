import os
import json
from unittest import TestCase
from unittest.mock import patch, MagicMock
from convert_gvf_to_vcf.gather_metadata import gather_metadata, add_file_metadata


class TestGatherMetadata(TestCase):
    def setUp(self):
        input_folder = os.path.dirname(__file__)
        self.config = os.path.join(input_folder, "input", "test.config")
        self.study_accession = "estd22"
        self.assembly = os.path.join(input_folder, "input", "zebrafish.fa")
        self.assembly_report = os.path.join(input_folder, "input", "assembly_report.txt")

        self.json_file = os.path.join(input_folder, "output", "a_preconverted.json")
        self.json_file_preconverted = os.path.join(input_folder, "output", "a_preconverted_preconverted.json")
        self.vcf_output = os.path.join(input_folder, "output", "a.vcf")
        self.expected_json_output = os.path.join(input_folder, "output", "a.json")

    @patch('convert_gvf_to_vcf.gather_metadata.DGVaMetadataRetriever')
    def test_gather_metadata(self, MockRetriever):
        mock_instance = MockRetriever.return_value
        mock_instance.__enter__.return_value = mock_instance
        gather_metadata(
            self.config,
            self.json_file,
            self.study_accession,
            self.assembly,
            self.assembly_report
        )

        with open(self.json_file_preconverted, 'r') as f_out:
            metadata = json.load(f_out)

        self.assertIn(metadata["files"][0].get("fileSize"), [None,""])
        self.assertIn(metadata["files"][0].get("md5"), [None, ""])

    def test_add_file_metadata(self):
        add_file_metadata(self.config, self.json_file, self.vcf_output)

        with open(self.json_file, 'r') as f:
            metadata = json.load(f)
        with open(self.expected_json_output, 'r') as f_out:
            expected_metadata = json.load(f_out)
        self.assertEqual(metadata["files"][0]["fileName"], expected_metadata["files"][0]["fileName"])
        self.assertEqual(metadata["files"][0]["fileSize"], expected_metadata["files"][0]["fileSize"])
        self.assertEqual(metadata["files"][0]["md5"], expected_metadata["files"][0]["md5"])


