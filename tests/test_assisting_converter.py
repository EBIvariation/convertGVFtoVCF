#TODO: 1 tests
import os
import unittest

from convert_gvf_to_vcf.assistingconverter import get_gvf_attributes, generate_custom_structured_meta_line, \
    convert_gvf_attributes_to_vcf_values
from convert_gvf_to_vcf.utils import read_yaml, generate_symbolic_allele_dict
from convert_gvf_to_vcf.lookup import Lookup

class TestAssistingConverter(unittest.TestCase):
    def setUp(self):
        input_folder = os.path.dirname(__file__)
        self.input_file = os.path.join(input_folder, "input", "zebrafish.gvf")
        self.input_folder_parent = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'convert_gvf_to_vcf'))
        # the inputs below are INFO attribute files
        self.etc_folder =  os.path.join(self.input_folder_parent, "etc")
        # self.mapping_attribute_dict = read_yaml(
        #     os.path.join(self.etc_folder, 'attribute_mapper.yaml'))  # formerly attributes_mapper and INFOattributes
        # self.symbolic_allele_dictionary = generate_symbolic_allele_dict(self.mapping_attribute_dict)
        self.assembly = os.path.join(input_folder, "input", "zebrafish.fa")
        self.etc_folder =  os.path.join(self.input_folder_parent, "etc")
        self.output_file = os.path.join(input_folder, "input", "a.vcf")
        self.reference_lookup = Lookup(self.assembly)

    def test_generate_custom_structured_meta_line(self):
        field = "INFO"
        idkey = "END"
        number = "1"
        data_type = "Integer"
        description = "Testing the end position"
        # Testing None
        optional_data_none = None
        output_for_none = generate_custom_structured_meta_line(
            field,
            idkey,
            number,
            data_type,
            description,
            optional_data_none
        )
        expected_output_for_none = "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Testing the end position\">"
        assert output_for_none == expected_output_for_none
        # Testing optional data
        optional_data = {"test1": "A"}
        output_for_optional = generate_custom_structured_meta_line(
            field,
            idkey,
            number,
            data_type,
            description,
            optional_data
        )
        expected_output_for_optional = "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Testing the end position\",test1=\"A\">"
        assert output_for_optional == expected_output_for_optional

    def test_get_gvf_attributes(self):
        col9 = "ID=1;Name=nssv1412199;Alias=CNV28955;parent=nsv811094;Start_range=.,1;End_range=2,.;sample_name=Wilds2-3"
        gvf_attribute_dictionary = get_gvf_attributes(col9)
        expected_dictionary = {'ID': '1', 'Name': 'nssv1412199', 'Alias': 'CNV28955', 'parent': 'nsv811094', 'Start_range': ['.', '1'], 'End_range': ['2', '.'], 'sample_name': 'Wilds2-3'}
        assert gvf_attribute_dictionary == expected_dictionary

    def test_convert_gvf_attributes_to_vcf_values(self):
        #TODO: complete

        # col9 = "ID=1;Name=nssv1412199;Alias=CNV28955;parent=nsv811094;Start_range=.,1;End_range=2,.;sample_name=Wilds2-3;Genotype=0:1"
        # # self.reference_lookup.mapping_attribute_dict
        # field_lines_dictionary = {'ALT': [], 'INFO': [], 'FILTER': [], 'FORMAT': []}
        # all_possible_lines_dictionary = {'ALT': [], 'INFO': [], 'FILTER': [], 'FORMAT': []}
        # gvf_attribute_dictionary, vcf_info_values, vcf_format_values = convert_gvf_attributes_to_vcf_values(
        #     col9,
        #     self.reference_lookup.mapping_attribute_dict,
        #     field_lines_dictionary,
        #     all_possible_lines_dictionary
        # )
        # testing gvf_attribute_dictionary
        expected_gvf_attribute_dictionary = {'ID': '1', 'Name': 'nssv1412199', 'Alias': 'CNV28955', 'parent': 'nsv811094',
                               'Start_range': ['.', '1'], 'End_range': ['2', '.'], 'sample_name': 'Wilds2-3'}
        # assert gvf_attribute_dictionary == expected_gvf_attribute_dictionary
        # testing vcf_info_values
        # expected_info_values = {'ID': '1', 'NAME': 'nssv1412199', 'ALIAS': 'CNV28955'}
        # assert vcf_info_values == expected_info_values
        # testing vcf_format_values
        # expected_format_values = {}
