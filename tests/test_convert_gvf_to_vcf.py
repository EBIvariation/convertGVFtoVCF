import os.path
import unittest

from convert_gvf_to_vcf.convertGVFtoVCF import generate_custom_unstructured_metainfomation_line, read_in_gvf_file, \
    read_dgva_info_attributes, read_gvf_info_attributes, gvf_features_to_vcf_objects, format_vcf_datalines, \
    generate_vcf_metainformation, write_to_vcf_file, generate_all_possible_standard_structured_info_lines, \
    generate_all_possible_standard_structured_alt_lines, generate_all_possible_standard_structured_filter_lines, \
    generate_all_possible_standard_structured_format_lines
from convert_gvf_to_vcf.convertGVFtoVCF import VcfLine, GvfFeatureline


class TestConvertGVFtoVCF(unittest.TestCase):
    def setUp(self):
        input_folder = os.path.dirname(__file__)
        self.input_file = os.path.join(input_folder, "input", "zebrafish.gvf")
        self.input_folder_parent = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'convert_gvf_to_vcf'))
        # the inputs below are INFO attribute files
        self.dgva_input_file = os.path.join(self.input_folder_parent, "etc","dgvaINFOattributes.tsv")
        self.gvf_input_file = os.path.join(self.input_folder_parent, "etc","gvfINFOattributes.tsv")
        self.output_file = os.path.join(input_folder, "input", "a.vcf")

    def test_read_in_gvf_file(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        assert len(gvf_pragmas) > 1
        assert len(gvf_non_essential) > 1
        assert len(gvf_lines_obj_list) > 1

    def test_read_dgva_info_attributes(self):
        dgva_attribute_dict = read_dgva_info_attributes(self.dgva_input_file)
        assert len(dgva_attribute_dict) > 1

    def test_read_gvf_info_attributes(self):
        gvf_attribute_dict = read_gvf_info_attributes(self.gvf_input_file)
        assert len(gvf_attribute_dict) > 1

    def test_gvf_features_to_vcf_objects(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_dgva_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_gvf_info_attributes(self.gvf_input_file)
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_ALT = []
        lines_standard_INFO = []
        lines_standard_FILTER = []
        lines_standard_FORMAT = []
        # Dictionary for all possible VCF meta-information lines
        all_possible_ALT_lines = {}
        all_possible_INFO_lines = {}  # dictionary, ID => INFO meta-information line for that particular ID
        all_possible_FILTER_lines = {}
        all_possible_FORMAT_lines = {}  # dictionary, ID => FORMAT meta-information line for that particular ID
        vcf_data_lines, list_of_vcf_objects = gvf_features_to_vcf_objects(gvf_lines_obj_list, dgva_attribute_dict,
                                                                          gvf_attribute_dict, lines_custom_structured,
                                                                          lines_standard_ALT, lines_standard_INFO,
                                                                          lines_standard_FILTER, lines_standard_FORMAT,
                                                                          all_possible_ALT_lines,
                                                                          all_possible_INFO_lines,
                                                                          all_possible_FILTER_lines,
                                                                          all_possible_FORMAT_lines
                                                                          )
        assert len(vcf_data_lines) > 1
        assert len(list_of_vcf_objects) > 1

    def test_get_ref(self):
        gvf_feature_line = "1	DGVa	copy_number_loss	776614	786127	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=.,776614;End_range=786127,.;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."
        f_list = gvf_feature_line.split("\t")
        line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6], f_list[7], f_list[8])
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_dgva_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_gvf_info_attributes(self.gvf_input_file)
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_ALT = []
        lines_standard_INFO = []
        lines_standard_FILTER = []
        lines_standard_FORMAT = []
        # Dictionary for all possible VCF meta-information lines
        all_possible_ALT_lines = {}
        all_possible_INFO_lines = {}  # dictionary, ID => INFO meta-information line for that particular ID
        all_possible_FILTER_lines = {}
        all_possible_FORMAT_lines = {}  # dictionary, ID => FORMAT meta-information line for that particular ID

        v = VcfLine(line_object,
                    dgva_attribute_dict,
                    gvf_attribute_dict,
                    lines_custom_structured,
                    lines_standard_ALT,
                    lines_standard_INFO,
                    lines_standard_FILTER,
                    lines_standard_FORMAT,
                    all_possible_ALT_lines,
                    all_possible_INFO_lines,
                    all_possible_FILTER_lines,
                    all_possible_FORMAT_lines)
        reference_allele = v.get_ref()
        assert reference_allele == "."

    def test_print_vcf_datalines(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_dgva_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_gvf_info_attributes(self.gvf_input_file)
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_ALT = []
        lines_standard_INFO = []
        lines_standard_FILTER = []
        lines_standard_FORMAT = []
        # Dictionary for all possible VCF meta-information lines
        all_possible_ALT_lines = {}
        all_possible_INFO_lines = {}  # dictionary, ID => INFO meta-information line for that particular ID
        all_possible_FILTER_lines = {}
        all_possible_FORMAT_lines = {}  # dictionary, ID => FORMAT meta-information line for that particular ID
        vcf_data_lines, list_of_vcf_objects = gvf_features_to_vcf_objects(gvf_lines_obj_list,
                                                                          dgva_attribute_dict,
                                                                          gvf_attribute_dict,
                                                                          lines_custom_structured,
                                                                          lines_standard_ALT,
                                                                          lines_standard_INFO,
                                                                          lines_standard_FILTER,
                                                                          lines_standard_FORMAT,
                                                                          all_possible_ALT_lines,
                                                                          all_possible_INFO_lines,
                                                                          all_possible_FILTER_lines,
                                                                          all_possible_FORMAT_lines
                                                                          )
        formatted_vcf_datalines = format_vcf_datalines(list_of_vcf_objects)
        assert len(formatted_vcf_datalines) > 1

    def test_generate_vcf_metainformation(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_dgva_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_gvf_info_attributes(self.gvf_input_file)
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_ALT = []
        lines_standard_INFO = []
        lines_standard_FILTER = []
        lines_standard_FORMAT = []
        # Dictionary for all possible VCF meta-information lines
        all_possible_ALT_lines = {}
        all_possible_INFO_lines = {}  # dictionary, ID => INFO meta-information line for that particular ID
        all_possible_FILTER_lines = {}
        all_possible_FORMAT_lines = {}  # dictionary, ID => FORMAT meta-information line for that particular ID
        vcf_data_lines, list_of_vcf_objects = gvf_features_to_vcf_objects(gvf_lines_obj_list,
                                                                          dgva_attribute_dict,
                                                                          gvf_attribute_dict,
                                                                          lines_custom_structured,
                                                                          lines_standard_ALT,
                                                                          lines_standard_INFO,
                                                                          lines_standard_FILTER,
                                                                          lines_standard_FORMAT,
                                                                          all_possible_ALT_lines,
                                                                          all_possible_INFO_lines,
                                                                          all_possible_FILTER_lines,
                                                                          all_possible_FORMAT_lines
                                                                          )
        formatted_vcf_datalines = format_vcf_datalines(list_of_vcf_objects)
        lines_custom_unstructured = ['##fileformat=VCFv4.4','##fileDate=20150715', '##source=DGVa','##source=DGVa', '##genome-build=NCBI GRCz10']
        unique_pragmas_to_add = generate_vcf_metainformation(lines_custom_unstructured, gvf_pragmas, list_of_vcf_objects)
        assert len(unique_pragmas_to_add) > 1

    def test_generate_custom_unstructured_metainfomation_line(self):
        lines_custom_unstructured = ['##fileformat=VCFv4.4', '##fileDate=20150715', '##source=DGVa', '##source=DGVa',
                                     '##genome-build=NCBI GRCz10']
        formatted_string = generate_custom_unstructured_metainfomation_line("test_string_key", "test_string_value", lines_custom_unstructured)
        assert formatted_string == "##test_string_key=test_string_value"

    def test_write_to_vcf_file(self):
        test_string = "##fileformat=VCFv4.4"
        write_to_vcf_file(self.output_file, test_string)
        #self.assertIn(test_string, self.output_file)
        output_file_size = os.path.getsize(self.output_file)
        assert output_file_size > 0

    def test_generate_all_possible_standard_structured_info_lines(self):
        all_possible_INFO_lines = generate_all_possible_standard_structured_info_lines()
        assert len(all_possible_INFO_lines) > 0

    def test_generate_all_possible_standard_structured_alt_lines(self):
        all_possible_ALT_lines = generate_all_possible_standard_structured_alt_lines()
        assert len(all_possible_ALT_lines) > 0

    #TODO: uncomment this test once the function generate_all_possible_standard_structured_filter_lines has been filled in

    # def test_generate_all_possible_standard_structured_filter_lines(self):
    #     all_possible_FILTER_lines = generate_all_possible_standard_structured_filter_lines()
    #     assert len(all_possible_FILTER_lines) > 0

    def test_generate_all_possible_standard_structured_format_lines(self):
        all_possible_FORMAT_lines = generate_all_possible_standard_structured_format_lines()
        assert len(all_possible_FORMAT_lines) > 0

if __name__ == '__main__':
    unittest.main()
