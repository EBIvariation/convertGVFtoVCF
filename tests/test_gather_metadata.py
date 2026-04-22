import os
import json
from unittest import TestCase

from convert_gvf_to_vcf.gather_metadata import add_file_metadata


class TestGatherMetadata(TestCase):
    def setUp(self):
        input_folder = os.path.dirname(__file__)
        self.config = os.path.join(input_folder, "input", "test.config")

        self.json_file = os.path.join(input_folder, "output", "a_pre_conversion.json")
        self.vcf_output = os.path.join(input_folder, "output", "a.vcf")
        self.expected_json_output = os.path.join(input_folder, "output", "a.json")

    def test_add_file_metadata(self):
        add_file_metadata(self.config, self.json_file, self.vcf_output)

        with open(self.json_file, 'r') as f:
            metadata = json.load(f)
        with open(self.expected_json_output, 'r') as f_out:
            expected_metadata = json.load(f_out)
        self.assertEqual(metadata["files"][0]["fileName"], expected_metadata["files"][0]["fileName"])
        self.assertEqual(metadata["files"][0]["fileSize"], expected_metadata["files"][0]["fileSize"])
        self.assertEqual(metadata["files"][0]["md5"], expected_metadata["files"][0]["md5"])


