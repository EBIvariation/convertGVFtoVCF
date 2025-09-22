import os.path
import unittest

from convert_gvf_to_vcf.convertGVFtoVCF import generate_custom_unstructured_metainfomation_line, read_in_gvf_file, \
    gvf_features_to_vcf_objects, format_vcf_datalines, \
    generate_vcf_metainformation, generate_all_possible_standard_structured_info_lines, \
    generate_all_possible_standard_structured_alt_lines, generate_all_possible_standard_structured_filter_lines, \
    generate_all_possible_standard_structured_format_lines, generate_vcf_header_line, populate_sample_formats, \
    format_sample_values, read_info_attributes
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
        self.assembly = os.path.join(input_folder, "input", "zebrafish.fa")

    def test_read_in_gvf_file(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        assert len(gvf_pragmas) > 1
        assert len(gvf_non_essential) > 1
        assert len(gvf_lines_obj_list) > 1

    def test_read_info_attributes(self):
        dgva_attribute_dict = read_info_attributes(self.dgva_input_file)
        assert len(dgva_attribute_dict) > 1
        gvf_attribute_dict = read_info_attributes(self.gvf_input_file)
        assert len(gvf_attribute_dict) > 1

    def test_gvf_features_to_vcf_objects(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_info_attributes(self.gvf_input_file)
        assembly_file = self.assembly
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_ALT = []
        lines_standard_INFO = []
        lines_standard_FILTER = []
        lines_standard_FORMAT = []
        # Dictionary for all possible VCF meta-information lines
        all_possible_INFO_lines = generate_all_possible_standard_structured_info_lines()
        all_possible_ALT_lines = generate_all_possible_standard_structured_alt_lines()
        all_possible_FILTER_lines = generate_all_possible_standard_structured_filter_lines()
        all_possible_FORMAT_lines = generate_all_possible_standard_structured_format_lines()

        vcf_data_lines, list_of_vcf_objects = gvf_features_to_vcf_objects(gvf_lines_obj_list, dgva_attribute_dict,
                                                                          gvf_attribute_dict,
                                                                          assembly_file,
                                                                          lines_custom_structured,
                                                                          lines_standard_ALT, lines_standard_INFO,
                                                                          lines_standard_FILTER, lines_standard_FORMAT,
                                                                          all_possible_ALT_lines,
                                                                          all_possible_INFO_lines,
                                                                          all_possible_FILTER_lines,
                                                                          all_possible_FORMAT_lines
                                                                          )
        assert len(vcf_data_lines) > 1
        assert len(list_of_vcf_objects) > 1

    def test_adjust_pos(self):
        gvf_feature_line = "chromosome1	DGVa	copy_number_loss	77	78	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=.,776614;End_range=786127,.;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."
        f_list = gvf_feature_line.split("\t")
        line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6], f_list[7], f_list[8])
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_info_attributes(self.gvf_input_file)
        assembly_file = self.assembly
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_ALT = []
        lines_standard_INFO = []
        lines_standard_FILTER = []
        lines_standard_FORMAT = []
        # Dictionary for all possible VCF meta-information lines
        all_possible_INFO_lines = generate_all_possible_standard_structured_info_lines()
        all_possible_ALT_lines = generate_all_possible_standard_structured_alt_lines()
        all_possible_FILTER_lines = generate_all_possible_standard_structured_filter_lines()
        all_possible_FORMAT_lines = generate_all_possible_standard_structured_format_lines()

        v = VcfLine(line_object,
                    dgva_attribute_dict,
                    gvf_attribute_dict,
                    assembly_file,
                    lines_custom_structured,
                    lines_standard_ALT,
                    lines_standard_INFO,
                    lines_standard_FILTER,
                    lines_standard_FORMAT,
                    all_possible_ALT_lines,
                    all_possible_INFO_lines,
                    all_possible_FILTER_lines,
                    all_possible_FORMAT_lines)
        value_to_change = -7
        new_pos_value = v.adjust_pos(value_to_change)
        assert new_pos_value > 0

    def test_add_padded_base(self):
        gvf_feature_line = "chromosome1	DGVa	copy_number_loss	77	78	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=.,776614;End_range=786127,.;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."
        f_list = gvf_feature_line.split("\t")
        line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6], f_list[7], f_list[8])
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_info_attributes(self.gvf_input_file)
        assembly_file = self.assembly
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_ALT = []
        lines_standard_INFO = []
        lines_standard_FILTER = []
        lines_standard_FORMAT = []
        # Dictionary for all possible VCF meta-information lines
        all_possible_INFO_lines = generate_all_possible_standard_structured_info_lines()
        all_possible_ALT_lines = generate_all_possible_standard_structured_alt_lines()
        all_possible_FILTER_lines = generate_all_possible_standard_structured_filter_lines()
        all_possible_FORMAT_lines = generate_all_possible_standard_structured_format_lines()

        v = VcfLine(line_object,
                    dgva_attribute_dict,
                    gvf_attribute_dict,
                    assembly_file,
                    lines_custom_structured,
                    lines_standard_ALT,
                    lines_standard_INFO,
                    lines_standard_FILTER,
                    lines_standard_FORMAT,
                    all_possible_ALT_lines,
                    all_possible_INFO_lines,
                    all_possible_FILTER_lines,
                    all_possible_FORMAT_lines)
        print("before", v.pos, v.ref, v.alt)
        (padded_base, pos, ref, alt) = v.add_padded_base(False)
        print(padded_base, pos, ref, alt)

    def test_build_iupac_ambiguity_code(self):
        gvf_feature_line = "chromosome1	DGVa	copy_number_loss	77	78	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=.,776614;End_range=786127,.;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."
        f_list = gvf_feature_line.split("\t")
        line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6],
                                     f_list[7], f_list[8])
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_info_attributes(self.gvf_input_file)
        assembly_file = self.assembly
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_ALT = []
        lines_standard_INFO = []
        lines_standard_FILTER = []
        lines_standard_FORMAT = []
        # Dictionary for all possible VCF meta-information lines
        all_possible_INFO_lines = generate_all_possible_standard_structured_info_lines()
        all_possible_ALT_lines = generate_all_possible_standard_structured_alt_lines()
        all_possible_FILTER_lines = generate_all_possible_standard_structured_filter_lines()
        all_possible_FORMAT_lines = generate_all_possible_standard_structured_format_lines()

        v = VcfLine(line_object,
                    dgva_attribute_dict,
                    gvf_attribute_dict,
                    assembly_file,
                    lines_custom_structured,
                    lines_standard_ALT,
                    lines_standard_INFO,
                    lines_standard_FILTER,
                    lines_standard_FORMAT,
                    all_possible_ALT_lines,
                    all_possible_INFO_lines,
                    all_possible_FILTER_lines,
                    all_possible_FORMAT_lines)
        my_ipuac_dictionary = v.build_iupac_ambiguity_code()
        print("my_ipuac_dictionary", my_ipuac_dictionary)

    def test_convert_iupac_ambiguity_code(self):
        gvf_feature_line = "chromosome1	DGVa	copy_number_loss	77	78	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=.,776614;End_range=786127,.;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."
        f_list = gvf_feature_line.split("\t")
        line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6],
                                     f_list[7], f_list[8])
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_info_attributes(self.gvf_input_file)
        assembly_file = self.assembly
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_ALT = []
        lines_standard_INFO = []
        lines_standard_FILTER = []
        lines_standard_FORMAT = []
        # Dictionary for all possible VCF meta-information lines
        all_possible_INFO_lines = generate_all_possible_standard_structured_info_lines()
        all_possible_ALT_lines = generate_all_possible_standard_structured_alt_lines()
        all_possible_FILTER_lines = generate_all_possible_standard_structured_filter_lines()
        all_possible_FORMAT_lines = generate_all_possible_standard_structured_format_lines()

        v = VcfLine(line_object,
                    dgva_attribute_dict,
                    gvf_attribute_dict,
                    assembly_file,
                    lines_custom_structured,
                    lines_standard_ALT,
                    lines_standard_INFO,
                    lines_standard_FILTER,
                    lines_standard_FORMAT,
                    all_possible_ALT_lines,
                    all_possible_INFO_lines,
                    all_possible_FILTER_lines,
                    all_possible_FORMAT_lines)

        my_ipuac_dictionary = v.build_iupac_ambiguity_code()
        converted_ref_allele = v.convert_iupac_ambiguity_code(my_ipuac_dictionary)
        assert converted_ref_allele not in ["R", "Y", "M", "K", "S", "D", "W", "H", "B", "V", "D", "N"]


    def test_get_ref(self):
        gvf_feature_line = "chromosome1	DGVa	copy_number_loss	77	78	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=.,776614;End_range=786127,.;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."
        f_list = gvf_feature_line.split("\t")
        line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6], f_list[7], f_list[8])
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_info_attributes(self.gvf_input_file)
        assembly_file = self.assembly
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_ALT = []
        lines_standard_INFO = []
        lines_standard_FILTER = []
        lines_standard_FORMAT = []
        # Dictionary for all possible VCF meta-information lines
        all_possible_INFO_lines = generate_all_possible_standard_structured_info_lines()
        all_possible_ALT_lines = generate_all_possible_standard_structured_alt_lines()
        all_possible_FILTER_lines = generate_all_possible_standard_structured_filter_lines()
        all_possible_FORMAT_lines = generate_all_possible_standard_structured_format_lines()

        v = VcfLine(line_object,
                    dgva_attribute_dict,
                    gvf_attribute_dict,
                    assembly_file,
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
        assert len(reference_allele) != 0

    def test_generate_vcf_metainformation(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_info_attributes(self.gvf_input_file)
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_ALT = []
        lines_standard_INFO = []
        lines_standard_FILTER = []
        lines_standard_FORMAT = []
        # Dictionary for all possible VCF meta-information lines
        all_possible_INFO_lines = generate_all_possible_standard_structured_info_lines()
        all_possible_ALT_lines = generate_all_possible_standard_structured_alt_lines()
        all_possible_FILTER_lines = generate_all_possible_standard_structured_filter_lines()
        all_possible_FORMAT_lines = generate_all_possible_standard_structured_format_lines()

        vcf_data_lines, list_of_vcf_objects = gvf_features_to_vcf_objects(gvf_lines_obj_list,
                                                                          dgva_attribute_dict,
                                                                          gvf_attribute_dict,
                                                                          self.assembly,
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

        lines_custom_unstructured = ['##fileformat=VCFv4.4','##fileDate=20150715', '##source=DGVa','##source=DGVa', '##genome-build=NCBI GRCz10']
        unique_pragmas_to_add, samples = generate_vcf_metainformation(lines_custom_unstructured, gvf_pragmas,
                                                                      gvf_non_essential, list_of_vcf_objects)
        assert len(unique_pragmas_to_add) > 1 and len(samples) > 1

    def test_generate_vcf_metainformation(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_info_attributes(self.gvf_input_file)
        assembly_file = self.assembly
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_ALT = []
        lines_standard_INFO = []
        lines_standard_FILTER = []
        lines_standard_FORMAT = []
        # Dictionary for all possible VCF meta-information lines
        all_possible_INFO_lines = generate_all_possible_standard_structured_info_lines()
        all_possible_ALT_lines = generate_all_possible_standard_structured_alt_lines()
        all_possible_FILTER_lines = generate_all_possible_standard_structured_filter_lines()
        all_possible_FORMAT_lines = generate_all_possible_standard_structured_format_lines()

        vcf_data_lines, list_of_vcf_objects = gvf_features_to_vcf_objects(gvf_lines_obj_list,
                                                                          dgva_attribute_dict,
                                                                          gvf_attribute_dict,
                                                                          assembly_file,
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

        lines_custom_unstructured = ['##fileformat=VCFv4.4','##fileDate=20150715', '##source=DGVa','##source=DGVa', '##genome-build=NCBI GRCz10']
        unique_pragmas_to_add, samples = generate_vcf_metainformation(lines_custom_unstructured, gvf_pragmas,
                                                                      gvf_non_essential, list_of_vcf_objects)
        assert len(unique_pragmas_to_add) > 1 and len(samples) > 1

    def test_generate_vcf_header_line(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_info_attributes(self.gvf_input_file)
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_ALT = []
        lines_standard_INFO = []
        lines_standard_FILTER = []
        lines_standard_FORMAT = []
        # Dictionary for all possible VCF meta-information lines
        all_possible_INFO_lines = generate_all_possible_standard_structured_info_lines()
        all_possible_ALT_lines = generate_all_possible_standard_structured_alt_lines()
        all_possible_FILTER_lines = generate_all_possible_standard_structured_filter_lines()
        all_possible_FORMAT_lines = generate_all_possible_standard_structured_format_lines()

        vcf_data_lines, list_of_vcf_objects = gvf_features_to_vcf_objects(gvf_lines_obj_list,
                                                                          dgva_attribute_dict,
                                                                          gvf_attribute_dict,
                                                                          self.assembly,
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

        lines_custom_unstructured = ['##fileformat=VCFv4.4','##fileDate=20150715', '##source=DGVa','##source=DGVa', '##genome-build=NCBI GRCz10']
        unique_pragmas_to_add, samples = generate_vcf_metainformation(lines_custom_unstructured, gvf_pragmas,
                                                                      gvf_non_essential, list_of_vcf_objects)
        header_fields = generate_vcf_header_line(samples)
        assert len(header_fields) > 1

    def test_populate_sample_formats(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_info_attributes(self.gvf_input_file)
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_ALT = []
        lines_standard_INFO = []
        lines_standard_FILTER = []
        lines_standard_FORMAT = []
        # Dictionary for all possible VCF meta-information lines
        all_possible_INFO_lines = generate_all_possible_standard_structured_info_lines()
        all_possible_ALT_lines = generate_all_possible_standard_structured_alt_lines()
        all_possible_FILTER_lines = generate_all_possible_standard_structured_filter_lines()
        all_possible_FORMAT_lines = generate_all_possible_standard_structured_format_lines()

        vcf_data_lines, list_of_vcf_objects = gvf_features_to_vcf_objects(gvf_lines_obj_list,
                                                                          dgva_attribute_dict,
                                                                          gvf_attribute_dict,
                                                                          self.assembly,
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
        lines_custom_unstructured = ['##fileformat=VCFv4.4','##fileDate=20150715', '##source=DGVa','##source=DGVa', '##genome-build=NCBI GRCz10']
        unique_pragmas_to_add, samples = generate_vcf_metainformation(lines_custom_unstructured, gvf_pragmas,
                                                                      gvf_non_essential, list_of_vcf_objects)
        sample_name_format_value = populate_sample_formats(samples)
        assert len(sample_name_format_value) > 1


    def test_format_sample_values(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_info_attributes(self.gvf_input_file)
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_ALT = []
        lines_standard_INFO = []
        lines_standard_FILTER = []
        lines_standard_FORMAT = []
        # Dictionary for all possible VCF meta-information lines
        all_possible_INFO_lines = generate_all_possible_standard_structured_info_lines()
        all_possible_ALT_lines = generate_all_possible_standard_structured_alt_lines()
        all_possible_FILTER_lines = generate_all_possible_standard_structured_filter_lines()
        all_possible_FORMAT_lines = generate_all_possible_standard_structured_format_lines()
        vcf_data_lines, list_of_vcf_objects = gvf_features_to_vcf_objects(gvf_lines_obj_list,
                                                                          dgva_attribute_dict,
                                                                          gvf_attribute_dict,
                                                                          self.assembly,
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
        lines_custom_unstructured = ['##fileformat=VCFv4.4','##fileDate=20150715', '##source=DGVa','##source=DGVa', '##genome-build=NCBI GRCz10']
        unique_pragmas_to_add, samples = generate_vcf_metainformation(lines_custom_unstructured, gvf_pragmas,
                                                                      gvf_non_essential, list_of_vcf_objects)
        sample_name_format_value = populate_sample_formats(samples)
        sample_format_values_string = format_sample_values(sample_name_format_value)
        assert isinstance(sample_format_values_string, str)

    def test_format_vcf_datalines(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_info_attributes(self.gvf_input_file)
        assembly_file = self.assembly
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_ALT = []
        lines_standard_INFO = []
        lines_standard_FILTER = []
        lines_standard_FORMAT = []
        # Dictionary for all possible VCF meta-information lines
        all_possible_INFO_lines = generate_all_possible_standard_structured_info_lines()
        all_possible_ALT_lines = generate_all_possible_standard_structured_alt_lines()
        all_possible_FILTER_lines = generate_all_possible_standard_structured_filter_lines()
        all_possible_FORMAT_lines = generate_all_possible_standard_structured_format_lines()

        vcf_data_lines, list_of_vcf_objects = gvf_features_to_vcf_objects(gvf_lines_obj_list,
                                                                          dgva_attribute_dict,
                                                                          gvf_attribute_dict,
                                                                          assembly_file,
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
        lines_custom_unstructured = ['##fileformat=VCFv4.4','##fileDate=20150715', '##source=DGVa','##source=DGVa', '##genome-build=NCBI GRCz10']
        unique_pragmas_to_add, samples = generate_vcf_metainformation(lines_custom_unstructured, gvf_pragmas,
                                                                      gvf_non_essential, list_of_vcf_objects)
        formatted_vcf_datalines = format_vcf_datalines(list_of_vcf_objects, samples)
        assert len(formatted_vcf_datalines) > 1

    def test_generate_custom_unstructured_metainfomation_line(self):
        lines_custom_unstructured = ['##fileformat=VCFv4.4', '##fileDate=20150715', '##source=DGVa', '##source=DGVa',
                                     '##genome-build=NCBI GRCz10']
        formatted_string = generate_custom_unstructured_metainfomation_line("test_string_key", "test_string_value", lines_custom_unstructured)
        assert formatted_string == "##test_string_key=test_string_value"

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
