import unittest
import os
from convert_gvf_to_vcf.conversionstatistics import FileStatistics
from convert_gvf_to_vcf.convertGVFtoVCF import convert_gvf_pragmas_for_vcf_header, generate_vcf_header_line
from convert_gvf_to_vcf.lookup import Lookup


class TestConversionStatistics(unittest.TestCase):
    def setUp(self):
        input_folder = os.path.dirname(__file__)
        self.gvf_file = os.path.join(input_folder, "input", "zebrafish.gvf")
        self.vcf_file = os.path.join(input_folder, "input", "a.vcf")
        self.assembly = os.path.join(input_folder, "input", "zebrafish.fa")
        self.reference_lookup = Lookup(self.assembly)
        gvf_pragma = ['##gff-version 3', '##gvf-version 1.06', '##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=7955', '##file-date 2015-07-15', '##genome-build NCBI GRCz10']
        gvf_pragma_comments = ['#Study_accession: nstd62', '#Study_type: Control Set', '#Display_name: Brown_et_al_2012', '#Publication: PMID=22203992;Journal=Proceedings of the National Academy of Sciences of the United States of America;Paper_title=Extensive genetic diversity and substructuring among zebrafish strains revealed through copy number variant analysis.;Publication_year=2012', '#Study: First_author=Kim Brown;Description=Comparative genomic hybridization analysis of 3 laboratory and one wild zebrafish populations for Copy Number Variants', '#Assembly_name: GRCz10', '#subject: subject_name=Wilds2-3', '#subject: subject_name=Zon9', '#subject: subject_name=JenMale7;subject_sex=Male', '#subject: subject_name=JenMale6;subject_sex=Male', '#sample: sample_name=JenMale6;subject_name=JenMale6', '#sample: sample_name=Wilds2-3;subject_name=Wilds2-3', '#sample: sample_name=Zon9;subject_name=Zon9', '#sample: sample_name=JenMale7;subject_name=JenMale7', '#testing_unknown_pragma']
        unique_pragma, unique_sample_name =convert_gvf_pragmas_for_vcf_header(gvf_pragma, gvf_pragma_comments, self.reference_lookup)
        self.samples = unique_sample_name
        self.gvf_header = unique_pragma

    def test_get_file_md5(self):
        statistics_summariser = FileStatistics(gvf_file_path=self.gvf_file, gvf_pragmas=self.gvf_header, samples=self.samples)
        gvf_md5 = statistics_summariser.get_file_md5(self.gvf_file)
        assert gvf_md5 == "8d30ebd587b6b7d0d5475a23239f0748"

    def test_find_version(self):
        statistics_summariser = FileStatistics(gvf_file_path=self.gvf_file, gvf_pragmas=self.gvf_header, samples=self.samples)
        gvf_version = statistics_summariser.find_version(self.gvf_header)
        assert gvf_version == "1.06"

    def test___str__(self):
        statistics_summariser = FileStatistics(self.gvf_file, gvf_pragmas=self.gvf_header, samples=self.samples)
        statistics_summary = statistics_summariser.__str__()
        print(statistics_summary)
        (gvf_line,
         file_name_string,
         file_path,
         file_extension,
         file_size_string,
         timestamp_string,
         version,
         file_md5_string,
         gvf_chrom,
         vcf_chrom,
         gvf_feature_line_count,
         vcf_data_line_count,
         vcf_number_of_merges,
         vcf_alt_alleles_count,
         vcf_info_counter,
         vcf_format_counter,
         vcf_sample_number_count,
         vcf_alt_missing,
         vcf_imprecise_variants,
         self.vcf_precise_variants,
         end_string,
         _) = statistics_summary.split("\n")
        assert gvf_line == "======GVF======"
        assert file_name_string == "File name = zebrafish.gvf"
        assert file_extension == "File extension = .gvf"
        assert version == "Version = 1.06"
        # the below should be empty/set to 0 as gets set in convertGVFtoVCF:convert
        assert gvf_chrom == "gvf chrom = Counter()"
        assert vcf_chrom == "vcf chrom = Counter()"
        assert gvf_feature_line_count == "gvf_feature_line_count = 0"
        assert vcf_data_line_count == "vcf_data_line_count = 0"
        assert vcf_number_of_merges == "vcf_number_of_merges = 0"
        assert vcf_alt_alleles_count == "vcf_alt_allele_count = Counter()"
        assert vcf_info_counter == "vcf_info_counter = Counter()"
        assert vcf_sample_number_count == "vcf_sample_count = Counter()"
        assert end_string == "====="
