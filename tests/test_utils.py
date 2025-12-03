import os
import unittest

from convert_gvf_to_vcf.convertGVFtoVCF import generate_vcf_header_structured_lines
from convert_gvf_to_vcf.gvffeature import GvfFeatureline
from convert_gvf_to_vcf.utils import read_yaml, read_pragma_mapper, generate_symbolic_allele_dict, read_in_gvf_file, \
    build_iupac_ambiguity_code
from convert_gvf_to_vcf.vcfline import VcfLine
from convert_gvf_to_vcf.lookup import Lookup

class TestUtils(unittest.TestCase):
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


    def test_read_yaml(self):
        test_yaml_dictionary = read_yaml(os.path.join(self.etc_folder, 'attribute_mapper.yaml'))
        assert len(test_yaml_dictionary) > 0

    def test_read_pragma_mapper(self):
        pragma_to_vcf_header = read_pragma_mapper(os.path.join(self.etc_folder, 'pragma_mapper.tsv'))
        assert len(pragma_to_vcf_header) > 0

    def test_read_mapping_dictionary(self):
        symbolic_allele_dictionary = generate_symbolic_allele_dict(self.reference_lookup.mapping_attribute_dict)
        assert len(symbolic_allele_dictionary) > 0

    def test_read_in_gvf_file(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        assert len(gvf_pragmas) > 1
        assert len(gvf_non_essential) > 1
        assert len(gvf_lines_obj_list) > 1

    def test_build_iupac_ambiguity_code(self):
        expected_dictionary_iupac ={
            'R': ['A', 'G'],
            'Y': ['C', 'T'],
            'M': ['A', 'C'],
            'K': ['G', 'T'],
            'S': ['C', 'G'],
            'D': ['A', 'G', 'T'],
            'W': ['A', 'T'],
            'H': ['A', 'C', 'T'],
            'B': ['C', 'G', 'T'],
            'V': ['A', 'C', 'G'],
            'N': ['A', 'C', 'G', 'T']
        }
        dictionary_iupac = build_iupac_ambiguity_code()
        assert dictionary_iupac == expected_dictionary_iupac