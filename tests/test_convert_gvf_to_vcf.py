import os.path
import unittest


from convert_gvf_to_vcf.lookup import Lookup
from convert_gvf_to_vcf.convertGVFtoVCF import generate_vcf_header_unstructured_line, read_in_gvf_file, \
    convert_gvf_features_to_vcf_objects, convert_gvf_pragmas_for_vcf_header, generate_vcf_header_line, compare_vcf_objects, \
    determine_merge_or_keep_vcf_objects, merge_vcf_objects, parse_pragma, get_pragma_name_and_value, get_pragma_tokens, \
    keep_vcf_objects, get_sample_name_from_pragma, get_unique_sample_names, convert_gvf_pragmas_to_vcf_header, \
    convert_gvf_pragma_comment_to_vcf_header


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
        self.ordered_list_of_samples = ['JenMale6', 'Wilds2-3', 'Zon9', 'JenMale7']


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
        gvf_pragma = "##file-date 2015-07-15"
        delimiter = " "
        name, value = parse_pragma(gvf_pragma, delimiter)
        assert name, value == "##file-date, 2015-07-15"
        # testing: pragma has only a name, no value, should print warning
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

    def test_get_sample_names_from_pragma(self):
        pragma_name = "#sample"
        pragma_value = "sample_name=JenMale6;subject_name=JenMale6"
        list_of_samples = get_sample_name_from_pragma(pragma_name, pragma_value)
        assert isinstance(list_of_samples, str)
        expected_value = 'JenMale6'
        assert list_of_samples == expected_value

    def test_get_unique_sample_names(self):
        duplicated_samples = ["JenMale6", "JenMale6", "JenMale7"]
        unique_samples = get_unique_sample_names(duplicated_samples)
        expected_list = ["JenMale6", "JenMale7"]
        assert isinstance(unique_samples, list)
        assert unique_samples == expected_list

    def test_convert_gvf_pragmas_to_vcf_header(self):
        list_of_gvf_pragmas_to_convert = [
            '##gff-version 3',
            '##gvf-version 1.06',
            '##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=7955',
            '##file-date 2015-07-15',
            '##genome-build NCBI GRCz10'
        ]
        list_of_gvf_pragmas = ['##file-date', '##gff-version', '##gvf-version', '##species', '##genome-build']
        list_of_converted_pragmas = convert_gvf_pragmas_to_vcf_header(
            list_of_gvf_pragmas_to_convert,
            list_of_gvf_pragmas,
            self.reference_lookup.pragma_to_vcf_map
                                  )
        assert isinstance(list_of_converted_pragmas, list)
        expected_list = ['##gff-version=3', '##gvf-version=1.06', '##species=http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=7955', '##fileDate=2015-07-15', '##genome-build=NCBIGRCz10']
        assert list_of_converted_pragmas == expected_list

    def test_convert_gvf_pragma_comment_to_vcf_header(self):
        gvf_pragma_comments_to_convert = ['#Study_accession: nstd62', '#Study_type: Control Set', '#Display_name: Brown_et_al_2012', '#Publication: PMID=22203992;Journal=Proceedings of the National Academy of Sciences of the United States of America;Paper_title=Extensive genetic diversity and substructuring among zebrafish strains revealed through copy number variant analysis.;Publication_year=2012', '#Study: First_author=Kim Brown;Description=Comparative genomic hybridization analysis of 3 laboratory and one wild zebrafish populations for Copy Number Variants', '#Assembly_name: GRCz10', '#subject: subject_name=Wilds2-3', '#subject: subject_name=Zon9', '#subject: subject_name=JenMale7;subject_sex=Male', '#subject: subject_name=JenMale6;subject_sex=Male', '#sample: sample_name=JenMale6;subject_name=JenMale6', '#sample: sample_name=Wilds2-3;subject_name=Wilds2-3', '#sample: sample_name=Zon9;subject_name=Zon9', '#sample: sample_name=JenMale7;subject_name=JenMale7', '#testing_unknown_pragma']
        list_of_gvf_pragma_comments = ['#sample', '#Study_accession', '#Study_type', '#Display_name', '#Publication#Study', '#Assembly_name',
         '#subject']
        list_of_converted_pragma_comments, sample_names = convert_gvf_pragma_comment_to_vcf_header(gvf_pragma_comments_to_convert,
                                                                                           list_of_gvf_pragma_comments,
                                                                                           self.reference_lookup.pragma_to_vcf_map)
        expected_list_of_converted_pragma_comments = ['##Study_accession=nstd62', '##Study_type=Control Set', '##Display_name=Brown_et_al_2012', '##PMID=22203992', '##Journal=Proceedings of the National Academy of Sciences of the United States of America', '##Paper_title=Extensive genetic diversity and substructuring among zebrafish strains revealed through copy number variant analysis.', '##Publication_year=2012', '##First_author=Kim Brown', '##Description=Comparative genomic hybridization analysis of 3 laboratory and one wild zebrafish populations for Copy Number Variants', '##Assembly_name=GRCz10', '##subject=subject_name=Wilds2-3', '##subject=subject_name=Zon9', '##subject=subject_name=JenMale7;subject_sex=Male', '##subject=subject_name=JenMale6;subject_sex=Male', '##sample=sample_name=JenMale6;subject_name=JenMale6', '##sample=sample_name=Wilds2-3;subject_name=Wilds2-3', '##sample=sample_name=Zon9;subject_name=Zon9', '##sample=sample_name=JenMale7;subject_name=JenMale7']
        expected_list_of_sample_names = ['JenMale6', 'Wilds2-3', 'Zon9', 'JenMale7']
        assert isinstance(list_of_converted_pragma_comments, list)
        assert isinstance(sample_names, list)
        assert list_of_converted_pragma_comments == expected_list_of_converted_pragma_comments
        assert sample_names == expected_list_of_sample_names

    def test_convert_gvf_pragmas_for_vcf_header(self):
        gvf_pragma = ['##gff-version 3', '##gvf-version 1.06', '##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=7955', '##file-date 2015-07-15', '##genome-build NCBI GRCz10']
        gvf_pragma_comments = ['#Study_accession: nstd62', '#Study_type: Control Set', '#Display_name: Brown_et_al_2012', '#Publication: PMID=22203992;Journal=Proceedings of the National Academy of Sciences of the United States of America;Paper_title=Extensive genetic diversity and substructuring among zebrafish strains revealed through copy number variant analysis.;Publication_year=2012', '#Study: First_author=Kim Brown;Description=Comparative genomic hybridization analysis of 3 laboratory and one wild zebrafish populations for Copy Number Variants', '#Assembly_name: GRCz10', '#subject: subject_name=Wilds2-3', '#subject: subject_name=Zon9', '#subject: subject_name=JenMale7;subject_sex=Male', '#subject: subject_name=JenMale6;subject_sex=Male', '#sample: sample_name=JenMale6;subject_name=JenMale6', '#sample: sample_name=Wilds2-3;subject_name=Wilds2-3', '#sample: sample_name=Zon9;subject_name=Zon9', '#sample: sample_name=JenMale7;subject_name=JenMale7', '#testing_unknown_pragma']
        unique_converted_pragmas, unique_sample_name =convert_gvf_pragmas_for_vcf_header(gvf_pragma, gvf_pragma_comments, self.reference_lookup)
        expected_unique_converted_pragmas = ['##fileformat=VCFv4.4', '##gff-version=3', '##gvf-version=1.06', '##species=http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=7955', '##fileDate=2015-07-15', '##genome-build=NCBIGRCz10', '##Study_accession=nstd62', '##Study_type=Control Set', '##Display_name=Brown_et_al_2012', '##PMID=22203992', '##Journal=Proceedings of the National Academy of Sciences of the United States of America', '##Paper_title=Extensive genetic diversity and substructuring among zebrafish strains revealed through copy number variant analysis.', '##Publication_year=2012', '##First_author=Kim Brown', '##Description=Comparative genomic hybridization analysis of 3 laboratory and one wild zebrafish populations for Copy Number Variants', '##Assembly_name=GRCz10', '##subject=subject_name=Wilds2-3', '##subject=subject_name=Zon9', '##subject=subject_name=JenMale7;subject_sex=Male', '##subject=subject_name=JenMale6;subject_sex=Male', '##sample=sample_name=JenMale6;subject_name=JenMale6', '##sample=sample_name=Wilds2-3;subject_name=Wilds2-3', '##sample=sample_name=Zon9;subject_name=Zon9', '##sample=sample_name=JenMale7;subject_name=JenMale7']
        expected_unique_sample_name = ['JenMale6', 'Wilds2-3', 'Zon9', 'JenMale7']
        assert unique_converted_pragmas == expected_unique_converted_pragmas
        assert unique_sample_name == expected_unique_sample_name

    def test_generate_vcf_metainfo(self):
        gvf_pragmas, gvf_pragma_comments, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        (
            unique_converted_pragmas,
            unique_sample_names
        ) = convert_gvf_pragmas_for_vcf_header(
            gvf_pragmas, gvf_pragma_comments, self.reference_lookup
        )
        assert unique_converted_pragmas == ['##fileformat=VCFv4.4', '##gff-version=3', '##gvf-version=1.06', '##species=http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=7955', '##fileDate=2015-07-15', '##genome-build=NCBIGRCz10', '##Study_accession=nstd62', '##Study_type=Control Set', '##Display_name=Brown_et_al_2012', '##PMID=22203992', '##Journal=Proceedings of the National Academy of Sciences of the United States of America', '##Paper_title=Extensive genetic diversity and substructuring among zebrafish strains revealed through copy number variant analysis.', '##Publication_year=2012', '##First_author=Kim Brown', '##Description=Comparative genomic hybridization analysis of 3 laboratory and one wild zebrafish populations for Copy Number Variants', '##Assembly_name=GRCz10', '##subject=subject_name=Wilds2-3', '##subject=subject_name=Zon9', '##subject=subject_name=JenMale7;subject_sex=Male', '##subject=subject_name=JenMale6;subject_sex=Male', '##sample=sample_name=JenMale6;subject_name=JenMale6', '##sample=sample_name=Wilds2-3;subject_name=Wilds2-3', '##sample=sample_name=Zon9;subject_name=Zon9', '##sample=sample_name=JenMale7;subject_name=JenMale7']
        assert unique_sample_names == ['JenMale6', 'Wilds2-3', 'Zon9', 'JenMale7']

    def test_generate_vcf_header_line(self):
        header_fields = generate_vcf_header_line(['JenMale6', 'Wilds2-3', 'Zon9', 'JenMale7'])
        assert header_fields == '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tJenMale6\tWilds2-3\tZon9\tJenMale7'

    def test_gvf_features_to_vcf_objects(self):
        gvf_pragmas, gvf_pragma_comments, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        # standard structured meta-information lines for this VCF file
        (
            header_lines_for_this_vcf,
            list_of_vcf_objects
        ) = convert_gvf_features_to_vcf_objects(gvf_lines_obj_list, self.reference_lookup, self.ordered_list_of_samples)
        assert len(gvf_pragmas) > 1
        assert len(gvf_pragma_comments) > 1
        assert len(gvf_lines_obj_list) > 1
        assert len(header_lines_for_this_vcf) > 1
        assert len(list_of_vcf_objects) > 1

    def test_compare_vcf_objects(self):
        gvf_pragmas, gvf_pragma_comments, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        header_standard_lines_dictionary, list_of_vcf_objects = convert_gvf_features_to_vcf_objects(
            gvf_lines_obj_list,  self.reference_lookup, self.ordered_list_of_samples
        )
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

    def test_merge_vcf_objects(self):
        gvf_pragmas, gvf_pragma_comments, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        list_of_samples = ['JenMale6', 'Wilds2-3', 'Zon9', 'JenMale7']
        (header_standard_lines_dictionary,
         list_of_vcf_objects) = convert_gvf_features_to_vcf_objects(
                gvf_lines_obj_list,
            self.reference_lookup,
            list_of_samples)
        # use lines 4 and 5 of gvf file
        previous = list_of_vcf_objects[3] # line 4
        current = list_of_vcf_objects[4] #line 5
        merged_object = merge_vcf_objects(previous, current, list_of_samples)
        assert merged_object.chrom == "chromosome1"
        assert merged_object.pos == 127
        assert merged_object.id == "13;14"
        assert merged_object.ref == "GTACG"
        assert merged_object.alt == "<DUP>"
        assert merged_object.qual == "."
        assert merged_object.filter == "."
        assert merged_object.info_dict == {'VARSEQ': '.', 'SVCID': 'CNV6230,CNV5711', 'END': '131', 'ALIAS': 'CNV6230,CNV5711', 'AC': '3', 'VARCALLSOID': 'SO:0001742', 'REMAP': '.69625,.85344', 'NAME': 'nssv1389474,nssv1388955', 'SVLEN': '4'}


    def test_keep_vcf_objects(self):
        gvf_pragmas, gvf_pragma_comments, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        header_standard_lines_dictionary, list_of_vcf_objects = convert_gvf_features_to_vcf_objects(
                gvf_lines_obj_list, self.reference_lookup, self.ordered_list_of_samples)
        list_of_samples = ['JenMale6', 'Wilds2-3', 'Zon9', 'JenMale7']
        previous_object = list_of_vcf_objects[1]
        kept_object = keep_vcf_objects(previous_object, list_of_samples)
        assert kept_object == previous_object


    def test_determine_merge_or_keep_vcf_objects(self):
        gvf_pragmas, gvf_pragma_comments, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        header_standard_lines_dictionary, list_of_vcf_objects = convert_gvf_features_to_vcf_objects(
            gvf_lines_obj_list,  self.reference_lookup, self.ordered_list_of_samples
        )
        list_of_samples = ['JenMale6', 'Wilds2-3', 'Zon9', 'JenMale7']
        flags_for_list_of_vcf_objects = compare_vcf_objects(list_of_vcf_objects)
        merged_or_kept_objects = determine_merge_or_keep_vcf_objects(list_of_vcf_objects, flags_for_list_of_vcf_objects, list_of_samples)
        assert len(merged_or_kept_objects) == 5 # 3 kept + 2 merged
        # check variant 13 and 14 have been merged
        assert merged_or_kept_objects[3].id == "13;14"
        assert merged_or_kept_objects[3].info_dict["NAME"] == "nssv1389474,nssv1388955"

if __name__ == '__main__':
    unittest.main()
