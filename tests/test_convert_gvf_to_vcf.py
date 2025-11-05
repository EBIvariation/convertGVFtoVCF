import os.path
import unittest

#from convert_gvf_to_vcf.utils import read_file
from convert_gvf_to_vcf.convertGVFtoVCF import generate_custom_unstructured_meta_line, read_in_gvf_file, \
    gvf_features_to_vcf_objects, format_vcf_datalines, \
    generate_vcf_metainfo, generate_vcf_header_structured_lines,  \
    generate_vcf_header_line,  \
    format_sample_values, read_yaml, read_pragma_mapper, generate_symbolic_allele_dict

from convert_gvf_to_vcf.vcfline import VcfLine
from convert_gvf_to_vcf.gvffeature import GvfFeatureline

class TestConvertGVFtoVCF(unittest.TestCase):
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

    def test_read_yaml(self):
        test_yaml_dictionary = read_yaml(os.path.join(self.etc_folder, 'attribute_mapper.yaml'))
        assert len(test_yaml_dictionary) > 0

    def test_read_pragma_mapper(self):
        pragma_to_vcf_header = read_pragma_mapper(os.path.join(self.etc_folder, 'pragma_mapper.tsv'))
        assert len(pragma_to_vcf_header) > 0

    def test_read_mapping_dictionary(self):
        symbolic_allele_dictionary = generate_symbolic_allele_dict(self.mapping_attribute_dict)
        assert len(symbolic_allele_dictionary) > 0

    def test_read_in_gvf_file(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        assert len(gvf_pragmas) > 1
        assert len(gvf_non_essential) > 1
        assert len(gvf_lines_obj_list) > 1

    def test_gvf_features_to_vcf_objects(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        assembly_file = self.assembly

        # standard structured meta-information lines for this VCF file
        header_lines_for_this_vcf, vcf_data_lines, list_of_vcf_objects = gvf_features_to_vcf_objects(gvf_lines_obj_list,
                                                                          assembly_file,self.mapping_attribute_dict, self.symbolic_allele_dictionary)
        assert len(gvf_pragmas) > 1
        assert len(gvf_non_essential) > 1
        assert len(gvf_lines_obj_list) > 1
        assert len(header_lines_for_this_vcf) > 1
        assert len(vcf_data_lines) > 1
        assert len(list_of_vcf_objects) > 1

    def test_add_padded_base(self):
        gvf_feature_line = "chromosome1	DGVa	copy_number_loss	77	78	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=.,776614;End_range=786127,.;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."
        f_list = gvf_feature_line.split("\t")
        line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6], f_list[7], f_list[8])

        mapping_attribute_dict = self.mapping_attribute_dict

        symbolic_allele_dictionary = self.symbolic_allele_dictionary
        assembly_file = self.assembly

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
        all_possible_info_lines = generate_vcf_header_structured_lines("INFO", mapping_attribute_dict)
        all_possible_alt_lines = generate_vcf_header_structured_lines("ALT", mapping_attribute_dict)
        all_possible_filter_lines = generate_vcf_header_structured_lines("FILTER", mapping_attribute_dict)
        all_possible_format_lines = generate_vcf_header_structured_lines("FORMAT", mapping_attribute_dict)

        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }

        v = VcfLine(line_object,
                    mapping_attribute_dict,
                    symbolic_allele_dictionary,
                    assembly_file,
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

        mapping_attribute_dict = self.mapping_attribute_dict
        symbolic_allele_dictionary = self.symbolic_allele_dictionary
        assembly_file = self.assembly

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
        all_possible_alt_lines = generate_vcf_header_structured_lines("ALT", self.mapping_attribute_dict)
        all_possible_info_lines = generate_vcf_header_structured_lines("INFO", self.mapping_attribute_dict)
        all_possible_filter_lines = generate_vcf_header_structured_lines("FILTER", self.mapping_attribute_dict)
        all_possible_format_lines = generate_vcf_header_structured_lines("FORMAT", self.mapping_attribute_dict)

        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }
        v = VcfLine(line_object,
                    mapping_attribute_dict,
                    symbolic_allele_dictionary,
                    assembly_file,
                    standard_lines_dictionary,
                    all_possible_lines_dictionary)

        my_ipuac_dictionary = v.build_iupac_ambiguity_code()
        assert len(my_ipuac_dictionary) > 0

    def test_convert_iupac_ambiguity_code(self):
        gvf_feature_line = "chromosome1	DGVa	copy_number_loss	77	78	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=.,776614;End_range=786127,.;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."
        f_list = gvf_feature_line.split("\t")
        line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6],
                                     f_list[7], f_list[8])

        mapping_attribute_dict = self.mapping_attribute_dict
        symbolic_allele_dictionary = self.symbolic_allele_dictionary
        assembly_file = self.assembly

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
        all_possible_alt_lines = generate_vcf_header_structured_lines("ALT", self.mapping_attribute_dict)
        all_possible_info_lines = generate_vcf_header_structured_lines("INFO", self.mapping_attribute_dict)
        all_possible_filter_lines = generate_vcf_header_structured_lines("FILTER", self.mapping_attribute_dict)
        all_possible_format_lines = generate_vcf_header_structured_lines("FORMAT", self.mapping_attribute_dict)

        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }
        v = VcfLine(line_object,
                    mapping_attribute_dict,
                    symbolic_allele_dictionary,
                    assembly_file,
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

        symbolic_allele_dictionary = self.symbolic_allele_dictionary
        assembly_file = self.assembly

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
        all_possible_alt_lines = generate_vcf_header_structured_lines("ALT", self.mapping_attribute_dict)
        all_possible_info_lines = generate_vcf_header_structured_lines("INFO", self.mapping_attribute_dict)
        all_possible_filter_lines = generate_vcf_header_structured_lines("FILTER", self.mapping_attribute_dict)
        all_possible_format_lines = generate_vcf_header_structured_lines("FORMAT", self.mapping_attribute_dict)

        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }
        v = VcfLine(line_object,
                    self.mapping_attribute_dict,
                    symbolic_allele_dictionary,
                    assembly_file,
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

        symbolic_allele_dictionary = self.symbolic_allele_dictionary
        assembly_file = self.assembly

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
        all_possible_info_lines = generate_vcf_header_structured_lines("INFO", self.mapping_attribute_dict)
        all_possible_alt_lines = generate_vcf_header_structured_lines("ALT", self.mapping_attribute_dict)
        all_possible_filter_lines = generate_vcf_header_structured_lines("FILTER", self.mapping_attribute_dict)
        all_possible_format_lines = generate_vcf_header_structured_lines("FORMAT", self.mapping_attribute_dict)

        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }
        v = VcfLine(line_object,
                    self.mapping_attribute_dict,
                    symbolic_allele_dictionary,
                    assembly_file,
                    standard_lines_dictionary,
                    all_possible_lines_dictionary)
        reference_allele = v.get_ref()
        assert len(reference_allele) != 0
        assert reference_allele == 'TA'

    def test_generate_symbolic_allele(self):
        gvf_feature_line = "chromosome1	DGVa	copy_number_loss	77	81	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=77,78;End_range=80,81;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."
        f_list = gvf_feature_line.split("\t")
        line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6], f_list[7], f_list[8])
        symbolic_allele_dictionary = self.symbolic_allele_dictionary
        assembly_file = self.assembly

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
        all_possible_alt_lines = generate_vcf_header_structured_lines("ALT", self.mapping_attribute_dict)
        all_possible_info_lines = generate_vcf_header_structured_lines("INFO", self.mapping_attribute_dict)
        all_possible_filter_lines = generate_vcf_header_structured_lines("FILTER", self.mapping_attribute_dict)
        all_possible_format_lines = generate_vcf_header_structured_lines("FORMAT", self.mapping_attribute_dict)

        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }
        v = VcfLine(line_object,
                    self.mapping_attribute_dict,
                    symbolic_allele_dictionary,
                    assembly_file,
                    standard_lines_dictionary,
                    all_possible_lines_dictionary)
        (output_symbolic_allele, info_field, output_lines_standard_alt, output_lines_standard_info) = v.generate_symbolic_allele(standard_lines_dictionary, all_possible_lines_dictionary)
        assert output_symbolic_allele == '<DEL>'
        print(info_field)
        assert info_field == ['ID=1;NAME=nssv1412199;ALIAS=CNV28955;VARCALLSOID=SO:0001743;SVCID=CNV28955;REMAP=.98857;VARSEQ=.', 'END=81', 'SVLEN=4', 'IMPRECISE', 'CIPOS=0,1', 'CIEND=0,1', 'END=80', 'SVLEN=4', 'IMPRECISE', 'CIPOS=1,2', 'CIEND=1,2']

        assert output_lines_standard_alt == ['##ALT=<ID=DEL,Description="Deletion">', '##ALT=<ID=DEL,Description="Deletion">']
        print(output_lines_standard_info)
        assert output_lines_standard_info == ['##INFO=<ID=ID,Number=.,Type=String,Description="A unique identifier">', '##INFO=<ID=NAME,Number=.,Type=String,Description="Name">', '##INFO=<ID=ALIAS,Number=.,Type=String,Description="Secondary Name">', '##INFO=<ID=VARCALLSOID,Number=.,Type=String,Description="Variant call Sequence ontology ID">', '##INFO=<ID=SVCID,Number=.,Type=Integer,Description="submitter variant call ID">', '##INFO=<ID=REMAP,Number=.,Type=Float,Description="Remap score">', '##INFO=<ID=VARSEQ,Number=.,Type=String,Description="Alleles found in an individual (or group of individuals).">', '##INFO=<ID=END,Number=1,Type=Integer,Description="End position on CHROM (used with symbolic alleles; see below) or End position of the longest variant described in this record">', '##INFO=<ID=SVLEN,Number=A,Type=String,Description="Length of structural variant">', '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">', '##INFO=<ID=CIPOS,Number=.,Type=Integer,Description="Confidence interval around POS for symbolic structural variants">', '##INFO=<ID=CIEND,Number=.,Type=Integer,Description="Confidence interval around END for symbolic structural variants">', '##INFO=<ID=END,Number=1,Type=Integer,Description="End position on CHROM (used with symbolic alleles; see below) or End position of the longest variant described in this record">', '##INFO=<ID=SVLEN,Number=A,Type=String,Description="Length of structural variant">', '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">', '##INFO=<ID=CIPOS,Number=.,Type=Integer,Description="Confidence interval around POS for symbolic structural variants">', '##INFO=<ID=CIEND,Number=.,Type=Integer,Description="Confidence interval around END for symbolic structural variants">']



    def test_get_alt(self):
        gvf_feature_line = "chromosome1	DGVa	copy_number_loss	77	81	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=77,78;End_range=80,81;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."
        f_list = gvf_feature_line.split("\t")
        line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6], f_list[7], f_list[8])


        symbolic_allele_dictionary = self.symbolic_allele_dictionary
        assembly_file = self.assembly

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
        all_possible_alt_lines = generate_vcf_header_structured_lines("ALT", self.mapping_attribute_dict)
        all_possible_info_lines = generate_vcf_header_structured_lines("INFO", self.mapping_attribute_dict)
        all_possible_filter_lines = generate_vcf_header_structured_lines("FILTER", self.mapping_attribute_dict)
        all_possible_format_lines = generate_vcf_header_structured_lines("FORMAT", self.mapping_attribute_dict)

        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }
        v = VcfLine(line_object,
                    self.mapping_attribute_dict,
                    symbolic_allele_dictionary,
                    assembly_file,
                    standard_lines_dictionary,
                    all_possible_lines_dictionary)
        alt_allele = v.get_alt(standard_lines_dictionary, all_possible_lines_dictionary)
        assert alt_allele == '<DEL>'

    def test_generate_vcf_metainformation(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)

        header_standard_lines_dictionary, vcf_data_lines, list_of_vcf_objects = gvf_features_to_vcf_objects(gvf_lines_obj_list,
                                                                          self.assembly, self.mapping_attribute_dict, self.symbolic_allele_dictionary)
        print("standard lines", header_standard_lines_dictionary)
        (unique_pragmas_to_add, sample_names,
         unique_alt_lines_to_add, unique_info_lines_to_add,
         unique_filter_lines_to_add, unique_format_lines_to_add) = generate_vcf_metainfo(
            gvf_pragmas, gvf_non_essential, list_of_vcf_objects, header_standard_lines_dictionary
        )
        print(unique_pragmas_to_add)
        assert unique_pragmas_to_add == ['##fileformat=VCFv4.4', '##gff-version=3', '##source=DGVa', '##gvf-version=1.06',
                                         '##species=http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=7955', '##fileDate=2015-07-15',
                                         '##genome-build=NCBIGRCz10', '##Study_accession=nstd62', '##Study_type=Control Set',
                                         '##Display_name=Brown_et_al_2012', '##Publication_year=2012',
                                         '##Description=Comparative genomic hybridization analysis of 3 laboratory and one wild zebrafish populations for Copy Number Variants',
                                         '##Assembly_name=GRCz10', '##subject=subject_name=Wilds2-3', '##subject=subject_name=Zon9', '##subject=subject_name=JenMale7;subject_sex=Male',
                                         '##subject=subject_name=JenMale6;subject_sex=Male', '##sample=sample_name=JenMale6;subject_name=JenMale6', '##sample=sample_name=Wilds2-3;subject_name=Wilds2-3',
                                         '##sample=sample_name=Zon9;subject_name=Zon9', '##sample=sample_name=JenMale7;subject_name=JenMale7']


        assert unique_alt_lines_to_add == ['##ALT=<ID=DEL,Description="Deletion">', '##ALT=<ID=DUP,Description="Duplication">']
        assert unique_info_lines_to_add ==  ['##INFO=<ID=ID,Number=.,Type=String,Description="A unique identifier">', '##INFO=<ID=NAME,Number=.,Type=String,Description="Name">', '##INFO=<ID=ALIAS,Number=.,Type=String,Description="Secondary Name">', '##INFO=<ID=VARCALLSOID,Number=.,Type=String,Description="Variant call Sequence ontology ID">', '##INFO=<ID=SVCID,Number=.,Type=Integer,Description="submitter variant call ID">', '##INFO=<ID=REMAP,Number=.,Type=Float,Description="Remap score">', '##INFO=<ID=VARSEQ,Number=.,Type=String,Description="Alleles found in an individual (or group of individuals).">', '##INFO=<ID=END,Number=1,Type=Integer,Description="End position on CHROM (used with symbolic alleles; see below) or End position of the longest variant described in this record">', '##INFO=<ID=SVLEN,Number=A,Type=String,Description="Length of structural variant">', '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">', '##INFO=<ID=CIPOS,Number=.,Type=Integer,Description="Confidence interval around POS for symbolic structural variants">', '##INFO=<ID=CIEND,Number=.,Type=Integer,Description="Confidence interval around END for symbolic structural variants">', '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">', '##INFO=<ID=DBXREF,Number=.,Type=String,Description="A database cross-reference">', '##INFO=<ID=AD,Number=R,Type=Integer,Description="Total read depth for each allele">']

    def test_generate_vcf_header_line(self):
        header_fields = generate_vcf_header_line(['JenMale6', 'Wilds2-3', 'Zon9', 'JenMale7'])
        assert header_fields == '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tJenMale6\tWilds2-3\tZon9\tJenMale7'

    def test_format_sample_values(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        # standard structured meta-information lines for this VCF file
        header_standard_lines_dictionary, vcf_data_lines, list_of_vcf_objects = gvf_features_to_vcf_objects(gvf_lines_obj_list,
                                                                          self.assembly, self.mapping_attribute_dict, self.symbolic_allele_dictionary)
        (
            unique_pragmas_to_add, samples, unique_alt_lines_to_add, unique_info_lines_to_add,
            unique_filter_lines_to_add, unique_format_lines_to_add
        ) = generate_vcf_metainfo(gvf_pragmas, gvf_non_essential, list_of_vcf_objects,
                                  header_standard_lines_dictionary)
        for vcf_obj in list_of_vcf_objects:
            sample_name_dict_format_kv = vcf_obj.format_dict
            sample_format_values_string = format_sample_values(sample_name_dict_format_kv, samples)
            assert isinstance(sample_format_values_string, str)

    def test_format_vcf_datalines(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        header_standard_lines_dictionary, vcf_data_lines, list_of_vcf_objects = gvf_features_to_vcf_objects(gvf_lines_obj_list, self.assembly, self.mapping_attribute_dict, self.symbolic_allele_dictionary)
        (
            unique_pragmas_to_add, samples, unique_alt_lines_to_add, unique_info_lines_to_add,
            unique_filter_lines_to_add, unique_format_lines_to_add
         ) = generate_vcf_metainfo(gvf_pragmas, gvf_non_essential, list_of_vcf_objects, header_standard_lines_dictionary)
        formatted_vcf_datalines = format_vcf_datalines(list_of_vcf_objects, samples)
        print(formatted_vcf_datalines)
        assert formatted_vcf_datalines == ['chromosome1\t1\t1\tAC\t<DEL>\t.\t.\tID=1;NAME=nssv1412199;ALIAS=CNV28955;VARCALLSOID=SO:0001743;SVCID=CNV28955;REMAP=.98857;VARSEQ=.;END=1;SVLEN=1\t.\t.\t.\t.\t.', 'chromosome1\t76\t1\tTAA\t<DEL>\t.\t.\tID=1;NAME=nssv1412199;ALIAS=CNV28955;VARCALLSOID=SO:0001743;SVCID=CNV28955;REMAP=.98857;VARSEQ=.;END=78;SVLEN=1;IMPRECISE;CIPOS=776537,776837;CIEND=776537,776837\t.\t.\t.\t.\t.', 'chromosome1\t126\t12\tCGTACGGTACG\t<DEL>\t.\t.\tID=12;NAME=nssv1406143;ALIAS=CNV22899;VARCALLSOID=SO:0001743;SVCID=CNV22899;REMAP=.87402;VARSEQ=.;END=131;SVLEN=5\t.\t.\t.\t.\t.', 'chromosome1\t127\t13\tGTACGTACG\t<DUP>\t.\t.\tID=13;NAME=nssv1389474;ALIAS=CNV6230;VARCALLSOID=SO:0001742;SVCID=CNV6230;REMAP=.69625;VARSEQ=.;END=131;SVLEN=4\t.\t.\t.\t.\t.', 'chromosome1\t127\t14\tGTACGTACG\t<DUP>\t.\t.\tID=14;NAME=nssv1388955;ALIAS=CNV5711;VARCALLSOID=SO:0001742;SVCID=CNV5711;REMAP=.85344;VARSEQ=.;AC=3;END=131;SVLEN=4\t.\t.\t.\t.\t.', 'chromosome1\t127\t14\tGTT\t<DUP>\t.\t.\tID=14;NAME=nssv1388955;ALIAS=CNV5711;VARCALLSOID=SO:0001742;SVCID=CNV5711;REMAP=.85344;VARSEQ=.;AC=3;DBXREF=mydata;AD=3;END=128;SVLEN=1\tAD\t3\t.\t.\t.', 'chromosome1\t127\t14\tGTT\t<DUP>\t.\t.\tID=14;NAME=nssv1388955;ALIAS=CNV5711;VARCALLSOID=SO:0001742;SVCID=CNV5711;REMAP=.85344;VARSEQ=.;AC=3;DBXREF=mydata;AD=3;END=128;SVLEN=1\tAD:GT\t.\t.\t.\t3:0:1']

    def test_generate_custom_unstructured_metainfomation_line(self):
        formatted_string = generate_custom_unstructured_meta_line("test_string_key", "test_string_value")
        assert formatted_string == "##test_string_key=test_string_value"

if __name__ == '__main__':
    unittest.main()
