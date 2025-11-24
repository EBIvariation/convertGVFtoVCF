#TODO: 9 tests
import os
import unittest

from convert_gvf_to_vcf.convertGVFtoVCF import generate_vcf_header_structured_lines, convert_gvf_features_to_vcf_objects, \
    generate_vcf_header_metainfo
from convert_gvf_to_vcf.gvffeature import GvfFeatureline
from convert_gvf_to_vcf.utils import read_yaml, generate_symbolic_allele_dict, read_in_gvf_file
from convert_gvf_to_vcf.vcfline import VcfLine
from convert_gvf_to_vcf.lookup import Lookup


class TestVcfline(unittest.TestCase):
    def setUp(self):
        # Prepare Directories
        self.input_folder_parent = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'convert_gvf_to_vcf'))
        self.etc_folder =  os.path.join(self.input_folder_parent, "etc")
        input_folder = os.path.dirname(__file__)
        # Prepare Inputs
        self.input_file = os.path.join(input_folder, "input", "zebrafish.gvf")
        self.input_folder_parent = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'convert_gvf_to_vcf'))
        # Prepare Outputs
        self.output_file = os.path.join(input_folder, "input", "a.vcf")
        # Prepare References
        self.assembly = os.path.join(input_folder, "input", "zebrafish.fa")
        self.reference_lookup = Lookup(self.assembly)
        # Set up GVF line object
        gvf_feature_line = "chromosome1	DGVa	copy_number_loss	77	78	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=.,776614;End_range=786127,.;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."
        f_list = gvf_feature_line.split("\t")
        gvf_line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6], f_list[7], f_list[8])
        # Set up of data structures
        # Dictionary of standard structured meta-information lines for this VCF file
        lines_standard_alt = []
        lines_standard_info = []
        lines_standard_filter = []
        lines_standard_format = []
        self.standard_lines_dictionary = {
            "ALT": lines_standard_alt,
            "INFO": lines_standard_info,
            "FILTER": lines_standard_filter,
            "FORMAT": lines_standard_format,
        }

        # Dictionary for all possible VCF meta-information lines
        all_possible_info_lines = generate_vcf_header_structured_lines("INFO", self.reference_lookup.mapping_attribute_dict)
        all_possible_alt_lines = generate_vcf_header_structured_lines("ALT", self.reference_lookup.mapping_attribute_dict)
        all_possible_filter_lines = generate_vcf_header_structured_lines("FILTER", self.reference_lookup.mapping_attribute_dict)
        all_possible_format_lines = generate_vcf_header_structured_lines("FORMAT", self.reference_lookup.mapping_attribute_dict)
        self.all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }

        self.v = VcfLine(gvf_line_object,
                    self.standard_lines_dictionary,
                    self.all_possible_lines_dictionary,
                    self.reference_lookup)

    def test_add_padded_base(self):
        test_ref = "A"
        test_alt = "T"
        padded_base, pos, ref, alt = self.v.add_padded_base(test_ref, test_alt, True, self.assembly)
        assert padded_base is not None
        assert pos is not None
        assert ref is not None
        assert alt is not None



    def test_convert_iupac_ambiguity_code(self):
        ref_to_convert = "TAGD"
        converted_ref_allele = self.v.convert_iupac_ambiguity_code(self.reference_lookup.iupac_ambiguity_dictionary, ref_to_convert)
        assert converted_ref_allele not in ["R", "Y", "M", "K", "S", "D", "W", "H", "B", "V", "D", "N"]

    def test_check_ref(self):

        reference_allele_to_check = "TGCR"
        new_ref = self.v.check_ref(reference_allele_to_check, self.reference_lookup)
        iupac_code = ["R", "Y", "M", "K", "S", "D", "W", "H", "B", "V", "D", "N"]
        assert all(code not in new_ref for code in iupac_code)

    def test_get_ref(self):
        reference_allele = self.v.get_ref(self.reference_lookup)
        assert len(reference_allele) != 0
        assert reference_allele == 'TA'

    def test_generate_symbolic_allele(self):
        (output_symbolic_allele,
         info_field,
         output_lines_standard_alt,
         output_lines_standard_info) = self.v.generate_symbolic_allele(self.standard_lines_dictionary,
                                                                       self.all_possible_lines_dictionary,
                                                                       self.reference_lookup.symbolic_allele_dictionary)
        assert output_symbolic_allele == '<DEL>'
        assert info_field == ['ID=1;NAME=nssv1412199;ALIAS=CNV28955;VARCALLSOID=SO:0001743;SVCID=CNV28955;REMAP=.98857;VARSEQ=.', 'END=77', 'SVLEN=1', 'END=78', 'SVLEN=1']
        assert output_lines_standard_alt == ['##ALT=<ID=DEL,Description="Deletion">', '##ALT=<ID=DEL,Description="Deletion">']
        assert output_lines_standard_info == ['##INFO=<ID=ID,Number=.,Type=String,Description="A unique identifier">', '##INFO=<ID=NAME,Number=.,Type=String,Description="Name">', '##INFO=<ID=ALIAS,Number=.,Type=String,Description="Secondary Name">', '##INFO=<ID=VARCALLSOID,Number=.,Type=String,Description="Variant call Sequence ontology ID">', '##INFO=<ID=SVCID,Number=.,Type=Integer,Description="submitter variant call ID">', '##INFO=<ID=REMAP,Number=.,Type=Float,Description="Remap score">', '##INFO=<ID=VARSEQ,Number=.,Type=String,Description="Alleles found in an individual (or group of individuals).">', '##INFO=<ID=END,Number=1,Type=Integer,Description="End position on CHROM (used with symbolic alleles; see below) or End position of the longest variant described in this record">', '##INFO=<ID=SVLEN,Number=A,Type=String,Description="Length of structural variant">', '##INFO=<ID=END,Number=1,Type=Integer,Description="End position on CHROM (used with symbolic alleles; see below) or End position of the longest variant described in this record">', '##INFO=<ID=SVLEN,Number=A,Type=String,Description="Length of structural variant">']

    def test_get_alt(self):
        alt_allele = self.v.get_alt(self.standard_lines_dictionary, self.all_possible_lines_dictionary, self.reference_lookup)
        assert alt_allele == '<DEL>'

    def test__str__(self):
        pass

    def test_merge_and_add(self):
        # previous="1"
        # current ="2"
        # delimiter =";"
        # merged_string = merge_and_add(previous, current, delimiter)
        # assert len(merged_string) > 1
        pass

    def test_put_GT_format_key_first(self):
        pass

    def test_format_sample_values(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        # standard structured meta-information lines for this VCF file
        (
            header_standard_lines_dictionary, vcf_data_lines, list_of_vcf_objects
        ) = convert_gvf_features_to_vcf_objects(gvf_lines_obj_list, self.reference_lookup)
        (
            unique_pragmas_to_add, samples, unique_alt_lines_to_add, unique_info_lines_to_add,
            unique_filter_lines_to_add, unique_format_lines_to_add
        ) = generate_vcf_header_metainfo(gvf_pragmas, gvf_non_essential, list_of_vcf_objects,
                                         header_standard_lines_dictionary)
        for vcf_obj in list_of_vcf_objects:
            sample_name_dict_format_kv = vcf_obj.format_dict
            # sample_format_values_string = format_sample_values(sample_name_dict_format_kv, samples)
            sample_format_values_string = vcf_obj.combine_format_values_by_sample(sample_name_dict_format_kv, samples)
            assert isinstance(sample_format_values_string, str)
        number_of_tokens_should_have = len(samples)
        tokens= sample_format_values_string.split("\t")
        actual_number_of_tokens = len(tokens)
        assert actual_number_of_tokens == number_of_tokens_should_have, f"must have {number_of_tokens_should_have}"
        assert sample_format_values_string == ".:.\t.:.\t.:.\t0:1:3", "String must match expected value"

    def test_info_list_to_dict(self):
        pass

    def test_merge_info_dicts(self):
        pass

    def test_merge_info_string(self):
        pass

    def test_merge_format_keys(self):
        pass

    def test_merge(self):
        pass

    def test_keep(self):
        pass
