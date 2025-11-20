#TODO: 3 tests
import os
import unittest

from convert_gvf_to_vcf.utils import read_yaml, generate_symbolic_allele_dict


class TestAssistingConverter(unittest.TestCase):
    def setUp(self):
        input_folder = os.path.dirname(__file__)
        self.input_file = os.path.join(input_folder, "input", "zebrafish.gvf")
        self.input_folder_parent = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'convert_gvf_to_vcf'))
        # the inputs below are INFO attribute files
        self.etc_folder =  os.path.join(self.input_folder_parent, "etc")
        self.mapping_attribute_dict = read_yaml(
            os.path.join(self.etc_folder, 'attribute_mapper.yaml'))  # formerly attributes_mapper and INFOattributes
        self.etc_folder =  os.path.join(self.input_folder_parent, "etc")
        self.symbolic_allele_dictionary = generate_symbolic_allele_dict(self.mapping_attribute_dict)
        self.output_file = os.path.join(input_folder, "input", "a.vcf")
        self.assembly = os.path.join(input_folder, "input", "zebrafish.fa")

    def test_generate_custom_structured_meta_line(self):
        pass

    def test_get_gvf_attributes(self):
        pass

    def test_convert_gvf_attributes_to_vcf_values(self):
        pass