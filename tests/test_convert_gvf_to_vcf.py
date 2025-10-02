import os.path
import unittest

from convert_gvf_to_vcf.convertGVFtoVCF import generate_custom_unstructured_metainformation_line, read_in_gvf_file, \
    gvf_features_to_vcf_objects, format_vcf_datalines, \
    generate_vcf_metainformation, generate_all_possible_standard_structured_lines,  \
    generate_vcf_header_line, populate_sample_formats, \
    format_sample_values, read_info_attributes, read_sequence_ontology_symbolic_allele
from convert_gvf_to_vcf.convertGVFtoVCF import VcfLine, GvfFeatureline


class TestConvertGVFtoVCF(unittest.TestCase):
    def setUp(self):
        input_folder = os.path.dirname(__file__)
        self.input_file = os.path.join(input_folder, "input", "zebrafish.gvf")
        self.input_folder_parent = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'convert_gvf_to_vcf'))
        # the inputs below are INFO attribute files
        self.dgva_input_file = os.path.join(self.input_folder_parent, "etc","dgvaINFOattributes.tsv")
        self.gvf_input_file = os.path.join(self.input_folder_parent, "etc","gvfINFOattributes.tsv")
        self.symbolic_allele_file = os.path.join(self.input_folder_parent,"etc", 'svALTkeys.tsv')
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

    def test_read_sequence_ontology_symbolic_allele(self):
        symbolic_allele_dictionary = read_sequence_ontology_symbolic_allele(self.symbolic_allele_file)
        assert len(symbolic_allele_dictionary) > 1

    def test_gvf_features_to_vcf_objects(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_info_attributes(self.gvf_input_file)
        symbolic_allele_dictionary = read_sequence_ontology_symbolic_allele(self.symbolic_allele_file)
        assembly_file = self.assembly
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_alt = []
        lines_standard_info = []
        lines_standard_filter = []
        lines_standard_format = []
        # merging
        standard_lines_dictionary = {
            "ALT": lines_standard_alt,
            "INFO": lines_standard_info,
            "FILTER": lines_standard_filter,
            "FORMAT": lines_standard_format,
        }
        # Dictionary for all possible VCF meta-information lines
        all_possible_info_lines = generate_all_possible_standard_structured_lines("INFO")
        all_possible_alt_lines = generate_all_possible_standard_structured_lines("ALT")
        all_possible_filter_lines = generate_all_possible_standard_structured_lines("FILTER")
        all_possible_format_lines = generate_all_possible_standard_structured_lines("FORMAT")
        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }
        vcf_data_lines, list_of_vcf_objects = gvf_features_to_vcf_objects(gvf_lines_obj_list, dgva_attribute_dict,
                                                                          gvf_attribute_dict,
                                                                          symbolic_allele_dictionary,
                                                                          assembly_file,
                                                                          lines_custom_structured,
                                                                          standard_lines_dictionary,
                                                                          all_possible_lines_dictionary
                                                                          )
        assert len(vcf_data_lines) > 1
        assert len(list_of_vcf_objects) > 1

    def test_add_padded_base(self):
        gvf_feature_line = "chromosome1	DGVa	copy_number_loss	77	78	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=.,776614;End_range=786127,.;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."
        f_list = gvf_feature_line.split("\t")
        line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6], f_list[7], f_list[8])
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_info_attributes(self.gvf_input_file)
        symbolic_allele_dictionary = read_sequence_ontology_symbolic_allele(self.symbolic_allele_file)
        assembly_file = self.assembly
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_alt = []
        lines_standard_info = []
        lines_standard_filter = []
        lines_standard_format = []
        # merging
        standard_lines_dictionary = {
            "ALT": lines_standard_alt,
            "INFO": lines_standard_info,
            "FILTER": lines_standard_filter,
            "FORMAT": lines_standard_format,
        }

        # Dictionary for all possible VCF meta-information lines
        all_possible_info_lines = generate_all_possible_standard_structured_lines("INFO")
        all_possible_alt_lines = generate_all_possible_standard_structured_lines("ALT")
        all_possible_filter_lines = generate_all_possible_standard_structured_lines("FILTER")
        all_possible_format_lines = generate_all_possible_standard_structured_lines("FORMAT")

        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }

        v = VcfLine(line_object,
                    dgva_attribute_dict,
                    gvf_attribute_dict,
                    symbolic_allele_dictionary,
                    assembly_file,
                    lines_custom_structured,
                    standard_lines_dictionary,
                    all_possible_lines_dictionary)

        test_ref = "A"
        test_alt = "T"
        padded_base, pos, ref, alt = v.add_padded_base(test_ref, test_alt, True)
        assert padded_base is not None
        assert pos is not None
        assert ref is not None
        assert alt is not None

    def test_build_iupac_ambiguity_code(self):
        gvf_feature_line = "chromosome1	DGVa	copy_number_loss	77	78	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=.,776614;End_range=786127,.;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."
        f_list = gvf_feature_line.split("\t")
        line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6],
                                     f_list[7], f_list[8])
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_info_attributes(self.gvf_input_file)
        symbolic_allele_dictionary = read_sequence_ontology_symbolic_allele(self.symbolic_allele_file)
        assembly_file = self.assembly
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_alt = []
        lines_standard_info = []
        lines_standard_filter = []
        lines_standard_format = []
        # merging
        standard_lines_dictionary = {
            "ALT": lines_standard_alt,
            "INFO": lines_standard_info,
            "FILTER": lines_standard_filter,
            "FORMAT": lines_standard_format,
        }

        # Dictionary for all possible VCF meta-information lines
        all_possible_alt_lines = generate_all_possible_standard_structured_lines("ALT")
        all_possible_info_lines = generate_all_possible_standard_structured_lines("INFO")
        all_possible_filter_lines = generate_all_possible_standard_structured_lines("FILTER")
        all_possible_format_lines = generate_all_possible_standard_structured_lines("FORMAT")

        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }
        v = VcfLine(line_object,
                    dgva_attribute_dict,
                    gvf_attribute_dict,
                    symbolic_allele_dictionary,
                    assembly_file,
                    lines_custom_structured,
                    standard_lines_dictionary,
                    all_possible_lines_dictionary)

        my_ipuac_dictionary = v.build_iupac_ambiguity_code()
        assert len(my_ipuac_dictionary) > 0

    def test_convert_iupac_ambiguity_code(self):
        gvf_feature_line = "chromosome1	DGVa	copy_number_loss	77	78	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=.,776614;End_range=786127,.;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."
        f_list = gvf_feature_line.split("\t")
        line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6],
                                     f_list[7], f_list[8])
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_info_attributes(self.gvf_input_file)
        symbolic_allele_dictionary = read_sequence_ontology_symbolic_allele(self.symbolic_allele_file)
        assembly_file = self.assembly
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_alt = []
        lines_standard_info = []
        lines_standard_filter = []
        lines_standard_format = []
        # merging
        standard_lines_dictionary = {
            "ALT": lines_standard_alt,
            "INFO": lines_standard_info,
            "FILTER": lines_standard_filter,
            "FORMAT": lines_standard_format,
        }
        # Dictionary for all possible VCF meta-information lines
        all_possible_alt_lines = generate_all_possible_standard_structured_lines("ALT")
        all_possible_info_lines = generate_all_possible_standard_structured_lines("INFO")
        all_possible_filter_lines = generate_all_possible_standard_structured_lines("FILTER")
        all_possible_format_lines = generate_all_possible_standard_structured_lines("FORMAT")

        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }
        v = VcfLine(line_object,
                    dgva_attribute_dict,
                    gvf_attribute_dict,
                    symbolic_allele_dictionary,
                    assembly_file,
                    lines_custom_structured,
                    standard_lines_dictionary,
                    all_possible_lines_dictionary)


        my_ipuac_dictionary = v.build_iupac_ambiguity_code()
        ref_to_convert = "TAGD"
        converted_ref_allele = v.convert_iupac_ambiguity_code(my_ipuac_dictionary, ref_to_convert)
        assert converted_ref_allele not in ["R", "Y", "M", "K", "S", "D", "W", "H", "B", "V", "D", "N"]

    def test_check_ref(self):
        gvf_feature_line = "chromosome1	DGVa	copy_number_loss	77	78	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=.,776614;End_range=786127,.;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."
        f_list = gvf_feature_line.split("\t")
        line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6],
                                     f_list[7], f_list[8])
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_info_attributes(self.gvf_input_file)
        symbolic_allele_dictionary = read_sequence_ontology_symbolic_allele(self.symbolic_allele_file)
        assembly_file = self.assembly
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_alt = []
        lines_standard_info = []
        lines_standard_filter = []
        lines_standard_format = []
        # merging
        standard_lines_dictionary = {
            "ALT": lines_standard_alt,
            "INFO": lines_standard_info,
            "FILTER": lines_standard_filter,
            "FORMAT": lines_standard_format,
        }
        # Dictionary for all possible VCF meta-information lines
        all_possible_alt_lines = generate_all_possible_standard_structured_lines("ALT")
        all_possible_info_lines = generate_all_possible_standard_structured_lines("INFO")
        all_possible_filter_lines = generate_all_possible_standard_structured_lines("FILTER")
        all_possible_format_lines = generate_all_possible_standard_structured_lines("FORMAT")

        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }
        v = VcfLine(line_object,
                    dgva_attribute_dict,
                    gvf_attribute_dict,
                    symbolic_allele_dictionary,
                    assembly_file,
                    lines_custom_structured,
                    standard_lines_dictionary,
                    all_possible_lines_dictionary)

        reference_allele_to_check = "TGCR"
        new_ref = v.check_ref(reference_allele_to_check)
        iupac_code = ["R", "Y", "M", "K", "S", "D", "W", "H", "B", "V", "D", "N"]
        assert all(code not in new_ref for code in iupac_code)


    def test_get_ref(self):
        gvf_feature_line = "chromosome1	DGVa	copy_number_loss	77	78	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=.,776614;End_range=786127,.;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."
        f_list = gvf_feature_line.split("\t")
        line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6], f_list[7], f_list[8])
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_info_attributes(self.gvf_input_file)
        symbolic_allele_dictionary = read_sequence_ontology_symbolic_allele(self.symbolic_allele_file)
        assembly_file = self.assembly
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_alt = []
        lines_standard_info = []
        lines_standard_filter = []
        lines_standard_format = []
        # merging
        standard_lines_dictionary = {
            "ALT": lines_standard_alt,
            "INFO": lines_standard_info,
            "FILTER": lines_standard_filter,
            "FORMAT": lines_standard_format,
        }
        # Dictionary for all possible VCF meta-information lines
        all_possible_info_lines = generate_all_possible_standard_structured_lines("INFO")
        all_possible_alt_lines = generate_all_possible_standard_structured_lines("ALT")
        all_possible_filter_lines = generate_all_possible_standard_structured_lines("FILTER")
        all_possible_format_lines = generate_all_possible_standard_structured_lines("FORMAT")

        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }
        v = VcfLine(line_object,
                    dgva_attribute_dict,
                    gvf_attribute_dict,
                    symbolic_allele_dictionary,
                    assembly_file,
                    lines_custom_structured,
                    standard_lines_dictionary,
                    all_possible_lines_dictionary)
        reference_allele = v.get_ref()
        assert len(reference_allele) != 0

    def test_generate_symbolic_allele(self):
        gvf_feature_line = "chromosome1	DGVa	copy_number_loss	77	81	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=77,78;End_range=80,81;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."
        f_list = gvf_feature_line.split("\t")
        line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6], f_list[7], f_list[8])
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_info_attributes(self.gvf_input_file)
        symbolic_allele_dictionary = read_sequence_ontology_symbolic_allele(self.symbolic_allele_file)
        assembly_file = self.assembly
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_alt = []
        lines_standard_info = []
        lines_standard_filter = []
        lines_standard_format = []
        # merging
        standard_lines_dictionary = {
            "ALT": lines_standard_alt,
            "INFO": lines_standard_info,
            "FILTER": lines_standard_filter,
            "FORMAT": lines_standard_format,
        }
        # Dictionary for all possible VCF meta-information lines
        all_possible_alt_lines = generate_all_possible_standard_structured_lines("ALT")
        all_possible_info_lines = generate_all_possible_standard_structured_lines("INFO")
        all_possible_filter_lines = generate_all_possible_standard_structured_lines("FILTER")
        all_possible_format_lines = generate_all_possible_standard_structured_lines("FORMAT")

        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }
        v = VcfLine(line_object,
                    dgva_attribute_dict,
                    gvf_attribute_dict,
                    symbolic_allele_dictionary,
                    assembly_file,
                    lines_custom_structured,
                    standard_lines_dictionary,
                    all_possible_lines_dictionary)
        output_symbolic_allele, self.info, output_lines_standard_ALT, output_lines_standard_INFO = v.generate_symbolic_allele(standard_lines_dictionary, all_possible_lines_dictionary)
        assert len(output_symbolic_allele) > 1

    def test_get_alt(self):
        gvf_feature_line = "chromosome1	DGVa	copy_number_loss	77	81	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=77,78;End_range=80,81;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."
        f_list = gvf_feature_line.split("\t")
        line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6], f_list[7], f_list[8])
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_info_attributes(self.gvf_input_file)
        symbolic_allele_dictionary = read_sequence_ontology_symbolic_allele(self.symbolic_allele_file)
        assembly_file = self.assembly
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_alt = []
        lines_standard_info = []
        lines_standard_filter = []
        lines_standard_format = []
        # merging
        standard_lines_dictionary = {
            "ALT": lines_standard_alt,
            "INFO": lines_standard_info,
            "FILTER": lines_standard_filter,
            "FORMAT": lines_standard_format,
        }
        # Dictionary for all possible VCF meta-information lines
        all_possible_alt_lines = generate_all_possible_standard_structured_lines("ALT")
        all_possible_info_lines = generate_all_possible_standard_structured_lines("INFO")
        all_possible_filter_lines = generate_all_possible_standard_structured_lines("FILTER")
        all_possible_format_lines = generate_all_possible_standard_structured_lines("FORMAT")

        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }
        v = VcfLine(line_object,
                    dgva_attribute_dict,
                    gvf_attribute_dict,
                    symbolic_allele_dictionary,
                    assembly_file,
                    lines_custom_structured,
                    standard_lines_dictionary,
                    all_possible_lines_dictionary)
        alt_allele = v.get_alt(standard_lines_dictionary, all_possible_lines_dictionary)
        assert len(alt_allele) > 0



    def test_generate_vcf_metainformation(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_info_attributes(self.gvf_input_file)
        symbolic_allele_dictionary = read_sequence_ontology_symbolic_allele(self.symbolic_allele_file)
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_alt = []
        lines_standard_info = []
        lines_standard_filter = []
        lines_standard_format = []
        # merging
        standard_lines_dictionary = {
            "ALT": lines_standard_alt,
            "INFO": lines_standard_info,
            "FILTER": lines_standard_filter,
            "FORMAT": lines_standard_format,
        }
        # Dictionary for all possible VCF meta-information lines
        all_possible_alt_lines = generate_all_possible_standard_structured_lines("ALT")
        all_possible_info_lines = generate_all_possible_standard_structured_lines("INFO")
        all_possible_filter_lines = generate_all_possible_standard_structured_lines("FILTER")
        all_possible_format_lines = generate_all_possible_standard_structured_lines("FORMAT")
        # merging
        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }
        vcf_data_lines, list_of_vcf_objects = gvf_features_to_vcf_objects(gvf_lines_obj_list,
                                                                          dgva_attribute_dict,
                                                                          gvf_attribute_dict,
                                                                          symbolic_allele_dictionary,
                                                                          self.assembly,
                                                                          lines_custom_structured,
                                                                          standard_lines_dictionary,
                                                                          all_possible_lines_dictionary
                                                                          )

        lines_custom_unstructured = ['##fileformat=VCFv4.4','##fileDate=20150715', '##source=DGVa','##source=DGVa', '##genome-build=NCBI GRCz10']
        unique_pragmas_to_add, samples = generate_vcf_metainformation(lines_custom_unstructured, gvf_pragmas,
                                                                      gvf_non_essential, list_of_vcf_objects)
        assert len(unique_pragmas_to_add) > 1 and len(samples) > 1

    def test_generate_vcf_metainformation(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_info_attributes(self.gvf_input_file)
        symbolic_allele_dictionary = read_sequence_ontology_symbolic_allele(self.symbolic_allele_file)
        assembly_file = self.assembly
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_alt = []
        lines_standard_info = []
        lines_standard_filter = []
        lines_standard_format = []
        # merging
        standard_lines_dictionary = {
            "ALT": lines_standard_alt,
            "INFO": lines_standard_info,
            "FILTER": lines_standard_filter,
            "FORMAT": lines_standard_format,
        }
        # Dictionary for all possible VCF meta-information lines
        all_possible_alt_lines = generate_all_possible_standard_structured_lines("ALT")
        all_possible_info_lines = generate_all_possible_standard_structured_lines("INFO")
        all_possible_filter_lines = generate_all_possible_standard_structured_lines("FILTER")
        all_possible_format_lines = generate_all_possible_standard_structured_lines("FORMAT")
        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }
        vcf_data_lines, list_of_vcf_objects = gvf_features_to_vcf_objects(gvf_lines_obj_list,
                                                                          dgva_attribute_dict,
                                                                          gvf_attribute_dict,
                                                                          symbolic_allele_dictionary,
                                                                          assembly_file,
                                                                          lines_custom_structured,
                                                                          standard_lines_dictionary,
                                                                          all_possible_lines_dictionary
                                                                          )

        lines_custom_unstructured = ['##fileformat=VCFv4.4','##fileDate=20150715', '##source=DGVa','##source=DGVa', '##genome-build=NCBI GRCz10']
        unique_pragmas_to_add, samples, unique_alt_lines_to_add, unique_info_lines_to_add, unique_filter_lines_to_add, unique_format_lines_to_add = generate_vcf_metainformation(lines_custom_unstructured, gvf_pragmas, gvf_non_essential, list_of_vcf_objects, standard_lines_dictionary)
        assert len(unique_pragmas_to_add) > 1 and len(samples) > 1

    def test_generate_vcf_header_line(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_info_attributes(self.gvf_input_file)
        symbolic_allele_dictionary = read_sequence_ontology_symbolic_allele(self.symbolic_allele_file)
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_alt = []
        lines_standard_info = []
        lines_standard_filter = []
        lines_standard_format = []
        # merging
        standard_lines_dictionary = {
            "ALT": lines_standard_alt,
            "INFO": lines_standard_info,
            "FILTER": lines_standard_filter,
            "FORMAT": lines_standard_format,
        }
        # Dictionary for all possible VCF meta-information lines
        all_possible_alt_lines = generate_all_possible_standard_structured_lines("ALT")
        all_possible_info_lines = generate_all_possible_standard_structured_lines("INFO")
        all_possible_filter_lines = generate_all_possible_standard_structured_lines("FILTER")
        all_possible_format_lines = generate_all_possible_standard_structured_lines("FORMAT")
        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }
        vcf_data_lines, list_of_vcf_objects = gvf_features_to_vcf_objects(gvf_lines_obj_list,
                                                                          dgva_attribute_dict,
                                                                          gvf_attribute_dict,
                                                                          symbolic_allele_dictionary,
                                                                          self.assembly,
                                                                          lines_custom_structured,
                                                                          standard_lines_dictionary,
                                                                          all_possible_lines_dictionary
                                                                          )

        lines_custom_unstructured = ['##fileformat=VCFv4.4','##fileDate=20150715', '##source=DGVa','##source=DGVa', '##genome-build=NCBI GRCz10']
        unique_pragmas_to_add, samples, unique_alt_lines_to_add, unique_info_lines_to_add, unique_filter_lines_to_add, unique_format_lines_to_add = generate_vcf_metainformation(lines_custom_unstructured, gvf_pragmas, gvf_non_essential, list_of_vcf_objects, standard_lines_dictionary)
        header_fields = generate_vcf_header_line(samples)
        assert len(header_fields) > 1

    def test_populate_sample_formats(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_info_attributes(self.gvf_input_file)
        symbolic_allele_dictionary = read_sequence_ontology_symbolic_allele(self.symbolic_allele_file)
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_alt = []
        lines_standard_info = []
        lines_standard_filter = []
        lines_standard_format = []
        # merging
        standard_lines_dictionary = {
            "ALT": lines_standard_alt,
            "INFO": lines_standard_info,
            "FILTER": lines_standard_filter,
            "FORMAT": lines_standard_format,
        }
        # Dictionary for all possible VCF meta-information lines

        all_possible_alt_lines = generate_all_possible_standard_structured_lines("ALT")
        all_possible_info_lines = generate_all_possible_standard_structured_lines("INFO")
        all_possible_filter_lines = generate_all_possible_standard_structured_lines("FILTER")
        all_possible_format_lines = generate_all_possible_standard_structured_lines("FORMAT")

        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }
        vcf_data_lines, list_of_vcf_objects = gvf_features_to_vcf_objects(gvf_lines_obj_list,
                                                                          dgva_attribute_dict,
                                                                          gvf_attribute_dict,
                                                                          symbolic_allele_dictionary,
                                                                          self.assembly,
                                                                          lines_custom_structured,
                                                                          standard_lines_dictionary,
                                                                          all_possible_lines_dictionary
                                                                          )
        lines_custom_unstructured = ['##fileformat=VCFv4.4','##fileDate=20150715', '##source=DGVa','##source=DGVa', '##genome-build=NCBI GRCz10']
        unique_pragmas_to_add, samples, unique_alt_lines_to_add, unique_info_lines_to_add, unique_filter_lines_to_add, unique_format_lines_to_add = generate_vcf_metainformation(
            lines_custom_unstructured, gvf_pragmas, gvf_non_essential, list_of_vcf_objects, standard_lines_dictionary)
        sample_name_format_value = populate_sample_formats(samples)
        assert len(sample_name_format_value) > 1


    def test_format_sample_values(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_info_attributes(self.gvf_input_file)
        symbolic_allele_dictionary = read_sequence_ontology_symbolic_allele(self.symbolic_allele_file)
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_alt = []
        lines_standard_info = []
        lines_standard_filter = []
        lines_standard_format = []
        # merging
        standard_lines_dictionary = {
            "ALT": lines_standard_alt,
            "INFO": lines_standard_info,
            "FILTER": lines_standard_filter,
            "FORMAT": lines_standard_format,
        }

        # Dictionary for all possible VCF meta-information lines
        all_possible_alt_lines = generate_all_possible_standard_structured_lines("ALT")
        all_possible_info_lines = generate_all_possible_standard_structured_lines("INFO")
        all_possible_filter_lines = generate_all_possible_standard_structured_lines("FILTER")
        all_possible_format_lines = generate_all_possible_standard_structured_lines("FORMAT")
        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }
        vcf_data_lines, list_of_vcf_objects = gvf_features_to_vcf_objects(gvf_lines_obj_list,
                                                                          dgva_attribute_dict,
                                                                          gvf_attribute_dict,
                                                                          symbolic_allele_dictionary,
                                                                          self.assembly,
                                                                          lines_custom_structured,
                                                                          standard_lines_dictionary,
                                                                          all_possible_lines_dictionary)

        lines_custom_unstructured = ['##fileformat=VCFv4.4','##fileDate=20150715', '##source=DGVa','##source=DGVa', '##genome-build=NCBI GRCz10']
        unique_pragmas_to_add, samples, unique_alt_lines_to_add, unique_info_lines_to_add, unique_filter_lines_to_add, unique_format_lines_to_add = generate_vcf_metainformation(
            lines_custom_unstructured, gvf_pragmas, gvf_non_essential, list_of_vcf_objects, standard_lines_dictionary)
        sample_name_format_value = populate_sample_formats(samples)
        sample_format_values_string = format_sample_values(sample_name_format_value)
        assert isinstance(sample_format_values_string, str)

    def test_format_vcf_datalines(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        dgva_attribute_dict = read_info_attributes(self.dgva_input_file)
        gvf_attribute_dict = read_info_attributes(self.gvf_input_file)
        symbolic_allele_dictionary = read_sequence_ontology_symbolic_allele(self.symbolic_allele_file)
        assembly_file = self.assembly
        # custom meta-information lines for this VCF fileT
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_alt = []
        lines_standard_info = []
        lines_standard_filter = []
        lines_standard_format = []
        # merging
        standard_lines_dictionary = {
            "ALT": lines_standard_alt,
            "INFO": lines_standard_info,
            "FILTER": lines_standard_filter,
            "FORMAT": lines_standard_format,
        }
        # Dictionary for all possible VCF meta-information lines
        all_possible_alt_lines = generate_all_possible_standard_structured_lines("ALT")
        all_possible_info_lines = generate_all_possible_standard_structured_lines("INFO")
        all_possible_filter_lines = generate_all_possible_standard_structured_lines("FILTER")
        all_possible_format_lines = generate_all_possible_standard_structured_lines("FORMAT")
        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }
        vcf_data_lines, list_of_vcf_objects = gvf_features_to_vcf_objects(gvf_lines_obj_list,
                                                                          dgva_attribute_dict,
                                                                          gvf_attribute_dict,
                                                                          symbolic_allele_dictionary,
                                                                          assembly_file,
                                                                          lines_custom_structured,
                                                                          standard_lines_dictionary,
                                                                          all_possible_lines_dictionary)
        lines_custom_unstructured = ['##fileformat=VCFv4.4','##fileDate=20150715', '##source=DGVa','##source=DGVa', '##genome-build=NCBI GRCz10']
        unique_pragmas_to_add, samples, unique_alt_lines_to_add, unique_info_lines_to_add, unique_filter_lines_to_add, unique_format_lines_to_add = generate_vcf_metainformation(lines_custom_unstructured, gvf_pragmas, gvf_non_essential, list_of_vcf_objects, standard_lines_dictionary)
        formatted_vcf_datalines = format_vcf_datalines(list_of_vcf_objects, samples)
        assert len(formatted_vcf_datalines) > 1

    def test_generate_custom_unstructured_metainfomation_line(self):
        lines_custom_unstructured = ['##fileformat=VCFv4.4', '##fileDate=20150715', '##source=DGVa', '##source=DGVa',
                                     '##genome-build=NCBI GRCz10']
        formatted_string = generate_custom_unstructured_metainformation_line("test_string_key", "test_string_value", lines_custom_unstructured)
        assert formatted_string == "##test_string_key=test_string_value"

if __name__ == '__main__':
    unittest.main()
