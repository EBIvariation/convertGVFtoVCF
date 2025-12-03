#TODO: 1 test
import os
import unittest

from convert_gvf_to_vcf.logger import set_up_logging


class TestLogger(unittest.TestCase):
    def setUp(self):
        input_folder = os.path.dirname(__file__)
        self.input_file = os.path.join(input_folder, "input", "zebrafish.gvf")
        self.input_folder_parent = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'convert_gvf_to_vcf'))
        # the inputs below are INFO attribute files
        self.etc_folder =  os.path.join(self.input_folder_parent, "etc")
        self.output_file = os.path.join(input_folder, "input", "a.vcf")
        self.assembly = os.path.join(input_folder, "input", "zebrafish.fa")

    def test_set_up_logging(self):
        # none as input
        expected_log_path_for_none = os.path.normpath(os.path.join(self.input_folder_parent, '..', "tests/output/converted.log"))
        log_path_for_none = set_up_logging(log_path=None)
        assert os.path.exists(log_path_for_none), f"{log_path_for_none} does not exist"
        assert log_path_for_none == expected_log_path_for_none
        # input path
        input_path = os.path.normpath(os.path.join(self.input_folder_parent, '..', 'tests/input/test.log'))
        log_path_user_input = set_up_logging(log_path=input_path)
        assert log_path_user_input == input_path
