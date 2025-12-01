import os.path
import unittest


from convert_gvf_to_vcf.lookup import Lookup
from convert_gvf_to_vcf.convertGVFtoVCF import generate_vcf_header_unstructured_line, read_in_gvf_file, \
    convert_gvf_features_to_vcf_objects, generate_vcf_header_metainfo, generate_vcf_header_line, compare_vcf_objects, \
    determine_merge_or_keep_vcf_objects, merge_vcf_objects, parse_pragma, get_pragma_name_and_value, get_pragma_tokens, \
    keep_vcf_objects



class TestConvertGVFtoVCF(unittest.TestCase):
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
        # self.mapping_attribute_dict = read_yaml(
        #     os.path.join(self.etc_folder, 'attribute_mapper.yaml'))  # formerly attributes_mapper and INFOattributes
        # self.symbolic_allele_dictionary = generate_symbolic_allele_dict(self.mapping_attribute_dict)



    def test_generate_vcf_header_structured_lines(self):
        key_to_test = "fileformat"
        value_to_test= "VCFv4.4"
        actual_output = generate_vcf_header_unstructured_line(key_to_test, value_to_test)
        assert actual_output == "##fileformat=VCFv4.4"

    def test_generate_custom_unstructured_meta_line(self):
        formatted_string = generate_vcf_header_unstructured_line("test_string_key", "test_string_value")
        assert formatted_string == "##test_string_key=test_string_value"

    def test_parse_pragma(self):
        # testing: pragma has a name and value
        essential_pragma = "##file-date 2015-07-15"
        delimiter = " "
        name, value = parse_pragma(essential_pragma, delimiter)
        assert name, value == "##file-date, 2015-07-15"
        # testing: pragma has only name, no value, should print warning
        name_only_pragma = "##file-date"
        name, value = parse_pragma(name_only_pragma, delimiter)
        assert name, value == "##file-date, None"
        # testing: invalid pragmas
        invalid_pragma = None
        with self.assertRaises(AttributeError):
            parse_pragma(invalid_pragma, delimiter)

    def test_get_pragma_name_and_value(self):
        pragma_to_test = "##file-date 2015-07-15"
        delimiter = " "
        list_of_pragma = ["##file-date", "##gff-version", "##gvf-version", "##species", "##genome-build"]
        pragma_to_vcf_map = {'##file-date': 'fileDate', '##gff-version': 'gff-version', '##gvf-version': 'gvf-version', '##species': 'species', '##genome-build': 'genome-build', '#sample': 'sample', '#Study_accession': 'Study_accession', '#Study_type': 'Study_type', '#Display_name': 'Display_name', '#Publication': 'Publication', '#Study': 'Study', '#Assembly_name': 'Assembly_name', '#subject': 'subject'}
        vcf_header_key, pragma_name, pragma_value = get_pragma_name_and_value(pragma_to_test, delimiter, list_of_pragma, pragma_to_vcf_map)
        assert vcf_header_key == "fileDate"
        assert pragma_name == "##file-date"
        assert pragma_value == "2015-07-15"

    def test_get_pragma_tokens(self):
        pragma_value = "First_author=Kim Brown;Description=Comparative genomic hybridization analysis of 3 laboratory and one wild zebrafish populations for Copy Number Variants"
        pragma_tokens = get_pragma_tokens(pragma_value, ";", "=")
        assert len(pragma_tokens) == 2
        # Testing: expected
        assert pragma_tokens[0][0] == "First_author"
        assert pragma_tokens[0][1] == "Kim Brown"
        assert pragma_tokens[1][0] == "Description"
        assert pragma_tokens[1][1] == "Comparative genomic hybridization analysis of 3 laboratory and one wild zebrafish populations for Copy Number Variants"
        # Testing: not expected
        unexpected_pragma_tokens = [['A', '1'], ['B', '2']]
        with self.assertRaises(AssertionError):
            self.assertEqual(unexpected_pragma_tokens, pragma_tokens)

    def test_generate_vcf_metainfo(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        (
            header_standard_lines_dictionary,
            vcf_data_lines,
            list_of_vcf_objects
        ) = convert_gvf_features_to_vcf_objects(gvf_lines_obj_list, self.reference_lookup)
        print("standard lines", header_standard_lines_dictionary)
        (unique_pragmas_to_add, sample_names,
         unique_alt_lines_to_add, unique_info_lines_to_add,
         unique_filter_lines_to_add, unique_format_lines_to_add) = generate_vcf_header_metainfo(
            gvf_pragmas, gvf_non_essential, list_of_vcf_objects, header_standard_lines_dictionary
        )
        print(unique_pragmas_to_add)
        assert unique_pragmas_to_add == ['##fileformat=VCFv4.4', '##gff-version=3', '##source=DGVa', '##gvf-version=1.06',
                                         '##species=http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=7955',
                                         '##fileDate=2015-07-15', '##genome-build=NCBIGRCz10', '##Study_accession=nstd62',
                                         '##Study_type=Control Set', '##Display_name=Brown_et_al_2012', '##PMID=22203992',
                                         '##Journal=Proceedings of the National Academy of Sciences of the United States of America',
                                         '##Paper_title=Extensive genetic diversity and substructuring among zebrafish strains revealed through copy number variant analysis.',
                                         '##Publication_year=2012', '##First_author=Kim Brown',
                                         '##Description=Comparative genomic hybridization analysis of 3 laboratory and one wild zebrafish populations for Copy Number Variants',
                                         '##Assembly_name=GRCz10', '##subject=subject_name=Wilds2-3', '##subject=subject_name=Zon9',
                                         '##subject=subject_name=JenMale7;subject_sex=Male', '##subject=subject_name=JenMale6;subject_sex=Male',
                                         '##sample=sample_name=JenMale6;subject_name=JenMale6', '##sample=sample_name=Wilds2-3;subject_name=Wilds2-3',
                                         '##sample=sample_name=Zon9;subject_name=Zon9', '##sample=sample_name=JenMale7;subject_name=JenMale7']



        assert unique_alt_lines_to_add == ['##ALT=<ID=DEL,Description="Deletion">', '##ALT=<ID=DUP,Description="Duplication">']
        assert unique_info_lines_to_add ==  ['##INFO=<ID=ID,Number=.,Type=String,Description="A unique identifier">', '##INFO=<ID=NAME,Number=.,Type=String,Description="Name">', '##INFO=<ID=ALIAS,Number=.,Type=String,Description="Secondary Name">', '##INFO=<ID=VARCALLSOID,Number=.,Type=String,Description="Variant call Sequence ontology ID">', '##INFO=<ID=SVCID,Number=.,Type=Integer,Description="submitter variant call ID">', '##INFO=<ID=REMAP,Number=.,Type=Float,Description="Remap score">', '##INFO=<ID=VARSEQ,Number=.,Type=String,Description="Alleles found in an individual (or group of individuals).">', '##INFO=<ID=END,Number=1,Type=Integer,Description="End position on CHROM (used with symbolic alleles; see below) or End position of the longest variant described in this record">', '##INFO=<ID=SVLEN,Number=A,Type=String,Description="Length of structural variant">', '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">', '##INFO=<ID=CIPOS,Number=.,Type=Integer,Description="Confidence interval around POS for symbolic structural variants">', '##INFO=<ID=CIEND,Number=.,Type=Integer,Description="Confidence interval around END for symbolic structural variants">', '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">', '##INFO=<ID=DBXREF,Number=.,Type=String,Description="A database cross-reference">', '##INFO=<ID=AD,Number=R,Type=Integer,Description="Total read depth for each allele">']

    def test_generate_vcf_header_line(self):
        header_fields = generate_vcf_header_line(['JenMale6', 'Wilds2-3', 'Zon9', 'JenMale7'])
        assert header_fields == '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tJenMale6\tWilds2-3\tZon9\tJenMale7'

    def test_gvf_features_to_vcf_objects(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        assembly_file = self.assembly

        # standard structured meta-information lines for this VCF file
        (
            header_lines_for_this_vcf,
            vcf_data_lines,
            list_of_vcf_objects
        ) = convert_gvf_features_to_vcf_objects(gvf_lines_obj_list, self.reference_lookup)
        assert len(gvf_pragmas) > 1
        assert len(gvf_non_essential) > 1
        assert len(gvf_lines_obj_list) > 1
        assert len(header_lines_for_this_vcf) > 1
        assert len(vcf_data_lines) > 1
        assert len(list_of_vcf_objects) > 1

    def test_compare_vcf_objects(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        print(gvf_lines_obj_list)
        header_standard_lines_dictionary, vcf_data_lines, list_of_vcf_objects = convert_gvf_features_to_vcf_objects(gvf_lines_obj_list, self.reference_lookup)
        print(list_of_vcf_objects)
        # compare object, if equal, True, if not equal, False # (next function will make true = current and merge; false= previous)
        expected_flags_for_list_of_vcf_objects = [False, # line 1 vs 2
                                                  False, # line 2 vs 3
                                                  False, # line 3 vs 4
                                                  True,  # line 4 vs 5
                                                  False, # line 5 vs 6
                                                  True   # line 6 vs 7
                                                  ]
        actual_flags_for_list_of_vcf_objects = compare_vcf_objects(list_of_vcf_objects)
        assert actual_flags_for_list_of_vcf_objects == expected_flags_for_list_of_vcf_objects

    # def test_merge_vcf_objects(self):
    #     gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
    #     header_standard_lines_dictionary, vcf_data_lines, list_of_vcf_objects = convert_gvf_features_to_vcf_objects(
    #         gvf_lines_obj_list, self.reference_lookup)
    #     list_of_samples = ['JenMale6', 'Wilds2-3', 'Zon9', 'JenMale7']
    #     # # use lines 4 and 5 of gvf file
    #     previous = list_of_vcf_objects[3] # line 4
    #     current = list_of_vcf_objects[4] #line 5
    #     merged_object = merge_vcf_objects(previous, current, list_of_samples)
    #
        # to_check = '\t'.join(['chromosome1', '127', '13;14', 'GTACGTACG', '<DUP>', '.', '.',
        #                       'ALIAS=CNV6230,CNV5711;NAME=nssv1389474,nssv1388955;VARSEQ=.;REMAP=.69625,.85344;SVCID=CNV6230,CNV5711;VARCALLSOID=SO:0001742;AC=3;SVLEN=4;END=131',
        #                       '.', '.\t.\t.\t.'])
        # assert merged_object == to_check

    def test_keep_vcf_objects(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        header_standard_lines_dictionary, vcf_data_lines, list_of_vcf_objects = convert_gvf_features_to_vcf_objects(
                gvf_lines_obj_list, self.reference_lookup)
        list_of_samples = ['JenMale6', 'Wilds2-3', 'Zon9', 'JenMale7']
        previous_object = list_of_vcf_objects[1]
        keep_vcf_objects(previous_object, list_of_samples)


    def test_determine_merge_or_keep_vcf_objects(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        header_standard_lines_dictionary, vcf_data_lines, list_of_vcf_objects = convert_gvf_features_to_vcf_objects(gvf_lines_obj_list, self.reference_lookup)
        list_of_samples = ['JenMale6', 'Wilds2-3', 'Zon9', 'JenMale7']
        flags_for_list_of_vcf_objects = compare_vcf_objects(list_of_vcf_objects)
        merged_or_kept_objects = determine_merge_or_keep_vcf_objects(list_of_vcf_objects, flags_for_list_of_vcf_objects, list_of_samples)
        for j in flags_for_list_of_vcf_objects:
            print(j)
        for i in merged_or_kept_objects:
            print(i)
        assert len(merged_or_kept_objects) == 5 # 3 kept + 2 merged
        # check variant 13 and 14 have been merged
        assert merged_or_kept_objects[3].id == "13;14"
        assert merged_or_kept_objects[3].info_dict["NAME"] == "nssv1389474,nssv1388955"

if __name__ == '__main__':
    unittest.main()
