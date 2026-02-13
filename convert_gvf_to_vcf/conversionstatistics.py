"""
The purpose of this file is to calculate the (conversion) statistics of the GVF and VCF file.
"""
import os
from datetime import datetime
import hashlib
from collections import Counter

class FileStatistics:
    """
    The responsibility of this class is to accumulate counts and calculate statistics of a GVF to VCF file conversion.
    """
    def __init__(self, gvf_file_path, gvf_pragmas, samples):
        self.gvf_file_path = gvf_file_path
        # file name
        self.gvf_file_name = os.path.basename(self.gvf_file_path)
        self.gvf_extension = os.path.splitext(self.gvf_file_path)[1]
        # file metadata
        self.gvf_metadata_stats = os.stat(self.gvf_file_path)
        self.gvf_file_size = self.gvf_metadata_stats.st_size
        gvf_last_modified_timestamp_raw = self.gvf_metadata_stats.st_mtime
        self.gvf_last_modified_timestamp = datetime.fromtimestamp(gvf_last_modified_timestamp_raw)
        self.gvf_version = self.find_version(gvf_pragmas)
        self.gvf_md5 = self.get_file_md5(self.gvf_file_path)
        # biological samples
        self.sample_number = len(samples) #
        self.vcf_sample_number_count = Counter()
        #TODO: add missing samples in VCF - awaiting bug fix
        #TODO: number of times a sample is missing in the merged vcf line, count per sample
        # chromsomes
        self.gvf_chromosome_count = Counter()
        self.vcf_chromosome_count = Counter()
        # # variants
        self.gvf_feature_line_count = 0
        self.gvf_sv_count = Counter()
        self.vcf_data_line_count = 0 # used in property
        self.vcf_number_of_merges = 0
        self.vcf_alt_alleles_count = Counter() #used in property
        # # attribute mapping
        self.vcf_info_counter = Counter() # used in property
        self.vcf_format_counter = Counter()
        self.vcf_variant_region_SOID = Counter()
        self.vcf_variant_call_SOID = Counter()

    @property
    def vcf_alt_missing(self):
        return self.vcf_alt_alleles_count['.']

    @property
    def vcf_imprecise_variants(self):
        return self.vcf_info_counter["IMPRECISE"]

    @property
    def vcf_precise_variants(self):
        return self.vcf_data_line_count - self.vcf_imprecise_variants

    @staticmethod
    def get_file_md5(path_to_file):
        hash_md5 = hashlib.md5()
        with open(path_to_file, "rb") as f:
            # Read file in 4KB chunks until you get to an empty byte string
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        md5 = hash_md5.hexdigest()
        return md5

    @staticmethod
    def find_version(header_list):
        """ Finds the gvf version in the header.
        :params: header_list
        :return: gvf version or if not found None
        """
        for line in header_list:
            if "##gvf-version" in line:
                try:
                    return line.split('=')[1].strip()
                except IndexError:
                    return f"GVF version unknown: {line}"
        return None

    def __str__(self):
        summary_string = (f"======GVF======\n"
                          f"File name = {self.gvf_file_name}\n"
                          f"File path = {self.gvf_file_path}\n"
                          f"File extension = {self.gvf_extension}\n"
                          f"File size = {self.gvf_file_size} bytes\n"
                          f"Last mod Timestamp = {self.gvf_last_modified_timestamp}\n"
                          f"Version = {self.gvf_version}\n"
                          f"md5 string = {self.gvf_md5}\n"
                          f"gvf chrom = {self.gvf_chromosome_count}\n"
                          f"vcf chrom = {self.vcf_chromosome_count}\n"
                          f"gvf_feature_line_count = {self.gvf_feature_line_count}\n"
                          f"vcf_data_line_count = {self.vcf_data_line_count}\n"
                          f"vcf_number_of_merges = {self.vcf_number_of_merges}\n"
                          f"vcf_alt_allele_count = {self.vcf_alt_alleles_count}\n"
                          f"vcf_info_counter = {self.vcf_info_counter}\n"
                          f"vcf_format_counter = {self.vcf_format_counter}\n"
                          f"vcf_sample_count = {self.vcf_sample_number_count}\n"
                          f"vcf_alt_missing = {self.vcf_alt_missing}\n"
                          f"vcf_imprecise_variants = {self.vcf_imprecise_variants}\n"
                          f"vcf_precise_variants = {self.vcf_precise_variants}\n"
                          f"=====\n"
                          )
        return summary_string

    def print_report(self, file_to_write):
        key_width = 30
        text_report = f'''
        GVF
            {"File name":<{key_width}}\t{self.gvf_file_name:>}
            {"File path":<{key_width}}\t{self.gvf_file_path:>}
            {"File extension":<{key_width}}\t{self.gvf_extension:>}
            {"File size":<{key_width}}\t{self.gvf_file_size:>}
            {"Last modified":<{key_width}}\t{self.gvf_last_modified_timestamp}
            {"Version":<{key_width}}\t{self.gvf_version:>}
            {"md5":<{key_width}}\t{self.gvf_md5:>}
            
             {"COUNTS":^{key_width}}
             {"GVF Chromosome counts:":<}\n\t\t\t\t{self.gvf_chromosome_count}
             {"Number of GVF feature lines:":<{key_width}}\t{self.gvf_feature_line_count:>}
             {"Number of samples:":<{key_width}}\t{self.sample_number}
             {"Number of structural variants seen:":<{key_width}}\n\t\t\t\t{self.gvf_sv_count}
        VCF
            {"COUNTS":^{key_width}}
            {"VCF Chromosome counts:":<{key_width}}\n\t\t\t\t{self.vcf_chromosome_count}
            {"VCF ALT allele counts:":<{key_width}}\n\t\t\t\t{self.vcf_alt_alleles_count}
            {"Number of VCF data lines:":<{key_width}}\t{self.vcf_data_line_count:>}
            {"Number of VCF data line merges:":<{key_width}}\t{self.vcf_number_of_merges:>}
            {"Number of times each VCF sample has been seen:":<{key_width}}\n\t\t\t\t{self.vcf_sample_number_count}
            {"Number of times INFO keys seen"}\n\t\t\t\t{self.vcf_info_counter}
            {"VCF variant region Sequence Ontology ID counts"}\n\t\t\t\t{self.vcf_variant_region_SOID}
            {"VCF variant call Sequence Ontology ID counts"}\n\t\t\t\t{self.vcf_variant_call_SOID}
            
        '''
        #
        with open(file_to_write, "w") as stats_file:
            stats_file.write(text_report)
