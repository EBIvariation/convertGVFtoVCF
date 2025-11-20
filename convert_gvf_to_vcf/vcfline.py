"""
The purpose of this file is to populate for each field of a VCF line (and perform any modifications/calculations to achieve this)
"""


from Bio import SeqIO
from convert_gvf_to_vcf.assistingconverter import convert_gvf_attributes_to_vcf_values


def extract_reference_allele(fasta_file, chromosome_name, position, end):
    """ Extracts the reference allele from the assembly.
    :param fasta_file: FASTA file of the assembly
    :param chromosome_name: name of the sequence
    :param position: position
    :param end: end position
    :return: reference_allele: base found at this chromosome_name at this position within this fasta_file
    """
    # using .index for memory efficiency:  https://biopython.org/docs/1.76/api/Bio.SeqIO.html#input-multiple-records
    records_dictionary = SeqIO.index(fasta_file, "fasta")
    zero_indexed_position = position - 1  # minus one because zero indexed
    zero_indexed_end = end - 1
    reference_allele = ""
    for position in range(zero_indexed_position, zero_indexed_end):
        reference_allele = reference_allele + records_dictionary[chromosome_name].seq[position]
    records_dictionary.close()
    return reference_allele

class VcfLine:
    """
    This class is responsible for the storing and merging of the fields of a VCF dataline.

    A VCF dataline is defined in the VCF specification as:
        - containing information about a position in the genome
        - genotype information on samples for each position.
    """
    def __init__(self,
                 gvf_feature_line_object,
                 field_lines_dictionary,
                 all_possible_lines_dictionary, #TODO: place this in reference
                 reference_lookup
                 ):

        (self.vcf_value, # used to populate the VCF fields. This is a dict of non-converted GVF attribute keys and their values.
         self.info_string, # a formatted INFO string to form VCF line
         self.format_dict # a dict of format tag and values for each sample to form VCF line
         ) = convert_gvf_attributes_to_vcf_values(gvf_feature_line_object.attributes, reference_lookup.mapping_attribute_dict, field_lines_dictionary, all_possible_lines_dictionary)

        # These might form useful parts of INFO field in VCF lines (useful information from GVF)
        self.source = gvf_feature_line_object.source
        self.so_type = gvf_feature_line_object.feature_type #currently column 3 of gvf, but could be an attribute so perhapsVCF: INFO or FORMAT?
        self.end = int(gvf_feature_line_object.end)
        self.phase = gvf_feature_line_object.phase # this is always a placeholder '.'

        # VCF DATALINE
        self.chrom = gvf_feature_line_object.seqid
        self.pos = int(gvf_feature_line_object.start)
        self.id = self.vcf_value["ID"]  # attributes: ID
        self.length = self.end - self.pos
        self.qual = gvf_feature_line_object.score # see EVA-3879: this is always '.'
        self.filter = "." # this is always a placeholder '.'; perhaps could add s50.

        # INFO
        self.key = self.chrom + "_" + str(self.pos)
        self.info = []
        self.info.append(self.info_string)
        # calculated last
        self.ref = self.get_ref(reference_lookup)
        self.alt = self.get_alt(field_lines_dictionary, all_possible_lines_dictionary, reference_lookup)

        self.sample_name = self.vcf_value["sample_name"] # this should be each samples names format value # sample names needs to be populated in attributes
        # # higher priority
        if self.format_dict:
            list_of_format_keys = [format_key for format_value in self.format_dict.values() for format_key in format_value.keys()]
            self.format = ":".join(list_of_format_keys)
        else:
            self.format = "." #TODO: this is temporary, when the multiple VCF lines are merged this will be filled in
        # list of samples from the VCF header should be here # merging will affect this # order of the sample

    # Functions which are responsible for token generation/population for the VCF line
    def add_padded_base(self, ref, alt, placed_before : bool, assembly_file):
        """ Adds a padded base to the REF and ALT allele of a VCF line.
        :param ref: reference allele
        :param alt: alt allele
        :param placed_before: padded base is placed before ref or alt True or False
        :return: (padded_base, self.pos, self.ref, self.alt)
        """
        if placed_before:
            padded_base_pos = self.pos - 1
            self.pos = padded_base_pos
            padded_base = extract_reference_allele(assembly_file, self.chrom, self.pos, self.end)
            ref = padded_base + ref
            if alt == ".":
                alt = padded_base
            else:
                alt = padded_base + alt
        elif not placed_before:
            padded_base_pos = self.pos + 1
            new_end = self.end + 1
            padded_base = extract_reference_allele(assembly_file, self.chrom, padded_base_pos, new_end)
            ref = ref + padded_base
            if alt == ".":
                alt = padded_base
            else:
                alt = alt + padded_base
        else:
            print("WARNING: Variable placed_before unknown: " + str(placed_before))
            padded_base = None
        return padded_base, self.pos, ref, alt

    def convert_iupac_ambiguity_code(self, iupac_ambiguity_dictionary, ref_to_convert):
        """ If the REF allele of a VCF line contains an IUPAC ambiguity code, converts it.
        :param iupac_ambiguity_dictionary: dictionary of IUPAC ambiguity code and a list of values
        :param ref_to_convert: reference allele to be converted
        :return: self.ref
        """
        converted_ref = ""
        for base in ref_to_convert:
            if base in iupac_ambiguity_dictionary:
                iupac_value = min(iupac_ambiguity_dictionary[base])
                converted_base = iupac_value
            else:
                converted_base = base
            converted_ref = converted_ref + converted_base
        return converted_ref

    def check_ref(self, ref_allele_to_be_checked, reference_lookup):
        """ Checks whether a reference allele meets the requirements of the VCF specification.
        :param ref_allele_to_be_checked: reference allele to check
        :return: checked_reference_allele: reference allele that meets the requirements of the VCF specification"""
        if isinstance(ref_allele_to_be_checked, str):
            if not all(bases in ref_allele_to_be_checked for bases in ["A", "C", "G", "T", "N"]):
                # checked_reference_allele = self.convert_iupac_ambiguity_code(self.build_iupac_ambiguity_code(), ref_allele_to_be_checked)
                checked_reference_allele = self.convert_iupac_ambiguity_code(reference_lookup.iupac_ambiguity_dictionary, ref_allele_to_be_checked)
            else:
                checked_reference_allele = ref_allele_to_be_checked
        else:
            print("WARNING: Ref allele must be a string. Setting to missing value.", ref_allele_to_be_checked)
            checked_reference_allele = "."
        return checked_reference_allele

    def get_ref(self, reference_lookup):
        """ Gets the reference allele from attributes column or if not found, returns "."
        :return: reference allele
        """
        assembly_file = reference_lookup.assembly_file
        if "Reference_seq" in self.vcf_value.keys():
            reference_allele = self.vcf_value["Reference_seq"]
        else:
            if assembly_file:
                reference_allele = extract_reference_allele(assembly_file, self.chrom, self.pos, self.end)
            else:
                print("WARNING: No reference provided. Placeholder inserted for Reference allele.")
                reference_allele = "."
        if reference_allele != ".":
            reference_allele = self.check_ref(reference_allele, reference_lookup)
        return reference_allele

    def generate_symbolic_allele(self, field_lines_dictionary, all_possible_lines_dictionary, symbolic_allele_dictionary):
        """ Generates the symbolic allele and stores the corresponding metainformation lines.
        Also determines if variant is precise or imprecise.
        :param field_lines_dictionary: lines for ALT, INFO, etc.
        :param all_possible_lines_dictionary: all possible lines
        :return: symbolic_allele, self.info, lines_standard_ALT, lines_standard_INFO
        """
        symbolic_allele_id = symbolic_allele_dictionary[self.so_type][1]
        symbolic_allele = f'<{symbolic_allele_id}>'

        lines_standard_alt = field_lines_dictionary["ALT"]
        lines_standard_info = field_lines_dictionary["INFO"]
        all_possible_alt_lines = all_possible_lines_dictionary["ALT"]
        all_possible_info_lines = all_possible_lines_dictionary["INFO"]

        if symbolic_allele_id in all_possible_alt_lines:
            lines_standard_alt.append(all_possible_alt_lines[symbolic_allele_id])

        info_svlen = None
        if self.length:
            info_svlen = "SVLEN=" + str(self.length)

        start_range_lower_bound = self.vcf_value["Start_range"][0]
        start_range_upper_bound = self.vcf_value["Start_range"][1]
        end_range_lower_bound = self.vcf_value["End_range"][0]
        end_range_upper_bound = self.vcf_value["End_range"][1]

        # setting up fields to be inserted into INFO
        info_end = None
        info_imprecise = None
        info_cipos = None
        info_ciend = None

        if start_range_lower_bound == "." or start_range_upper_bound == "." or end_range_lower_bound == "." or end_range_upper_bound == ".":
            is_imprecise = False
            info_end = "END=" + str(self.pos + len(self.ref) - 1)
        else:
            is_imprecise = True
            info_imprecise = "IMPRECISE"

            cipos_lower_bound = int(start_range_lower_bound) - self.pos
            cipos_upper_bound = int(start_range_upper_bound) - self.pos
            info_cipos = "CIPOS=" + str(cipos_lower_bound) + "," + str(cipos_upper_bound)

            ciend_lower_bound = int(start_range_lower_bound) - self.pos
            ciend_upper_bound = int(start_range_upper_bound) - self.pos
            info_ciend = "CIEND=" + str(ciend_lower_bound) + "," + str(ciend_upper_bound)

            if symbolic_allele == "<INS>":
                info_end ="END=" + str( self.pos + len(self.ref) - 1 )
            elif symbolic_allele in {"<DEL>", "<DUP>", "<INV>", "<CNV>"}:
                info_end = "END=" + str(self.pos + self.length)
            elif symbolic_allele == "<*>":
                info_end = "END=" + str(self.pos + len(self.ref))
            else:
                print("Cannot identify symbolic allele")

        # for all variants (precise and imprecise)
        self.info.append(info_end)
        lines_standard_info.append(all_possible_info_lines["END"])
        self.info.append(info_svlen)
        lines_standard_info.append(all_possible_info_lines["SVLEN"])

        # for imprecise variants only
        if is_imprecise:
            self.info.append(info_imprecise)
            lines_standard_info.append(all_possible_info_lines["IMPRECISE"])
            self.info.append(info_cipos)
            lines_standard_info.append(all_possible_info_lines["CIPOS"])
            self.info.append(info_ciend)
            lines_standard_info.append(all_possible_info_lines["CIEND"])
        return symbolic_allele, self.info, lines_standard_alt, lines_standard_info

    def get_alt(self, field_lines_dictionary, all_possible_lines_dictionary, reference_lookup):
        """ Gets the ALT allele for the VCF file
        :param field_lines_dictionary: store INFO,ALT, FILTER, FORMAT lines
        :param all_possible_lines_dictionary: dictionary of all possible ALT, INFO, FORMAT, FILTER lines
        :return: symbolic_allele, self.info, lines_standard_ALT, lines_standard_INFO
        """
        if any(base in self.vcf_value["Variant_seq"] for base in ["A", "C", "G", "T", "N"]):
            alterative_allele = self.vcf_value["Variant_seq"]
        elif self.vcf_value["Variant_seq"] == '.':
            symbolic_allele, self.info, lines_standard_alt, lines_standard_info = self.generate_symbolic_allele(field_lines_dictionary, all_possible_lines_dictionary, reference_lookup.symbolic_allele_dictionary)
            if symbolic_allele is None:
                alterative_allele = "."
            elif (self.vcf_value["Variant_seq"] == "." or self.vcf_value["Variant_seq"] == "-") and symbolic_allele is not None:
                alterative_allele = symbolic_allele
                # add padded bases
                if self.pos == 1:
                    #print("pos, ref, alt",self.pos,self.ref, alterative_allele)
                    padded_base, self.pos, self.ref, self.alt = self.add_padded_base(self.ref, alterative_allele, False, reference_lookup.assembly_file)
                    self.ref = self.check_ref(self.ref, reference_lookup)
                else:
                    #print("pos, ref, alt", self.pos,self.ref, alterative_allele)
                    padded_base, self.pos, self.ref, self.alt = self.add_padded_base(self.ref, alterative_allele, True, reference_lookup.assembly_file)
                    self.ref = self.check_ref(self.ref, reference_lookup)
            else:
                alterative_allele = "."
                print("Cannot identify symbolic allele. Variant type is not supported.")
        else:
            alterative_allele = "."
            print("Could not determine the alternative allele.")
        return alterative_allele

    def __str__(self):
        # string_to_return = '\t'.join((self.chrom, self.pos, self.key, self.qual, self.filter, self.info, self.source, self.phase, self.end, self.so_type, self.sample_name, self.format))
        self.info_string = self.format_info_string()
        string_to_return = '\t'.join((self.chrom,
                self.pos,
                self.id,
                self.ref,
                self.alt,
                self.qual,
                self.filter,
                self.info_string,
                self.format,
                self.format_values_by_sample_string
                ))
        return string_to_return

    def __eq__(self, other_vcf_line):
        """ Compares equality of PARTS of the VcfLine objects.
        :param: other_vcf_line: another object to compare equality with
        """
        if isinstance(other_vcf_line, VcfLine):
            return (self.chrom == other_vcf_line.chrom) and (self.pos == other_vcf_line.pos) and (self.ref == other_vcf_line.ref)
        return False

    def merge_and_add(self, previous_element, current_element, delimiter):
        """ Merges fields of a VCF line. If field is the same, use current element. If different, merge with delimiter.
        :param: previous_element
        :param: current_element
        :param: delimiter
        :return: merged element
        """
        if previous_element == current_element:
            merged_element = current_element
        else:
            merged_element = delimiter.join((previous_element, current_element))
        return merged_element
    # functions responsible for FORMAT are below
    def order_format_keys(self, set_of_format_keys):
        """Stores the FORMAT keys of the VCF line in the correct order by anchoring GT as the first key.
        :param: set_of_format_keys: format keys in a set
        :return: anchored_list_of_keys: list of ordered keys
        """
        anchored_list_of_format_keys = []
        if 'GT' in set_of_format_keys:
            anchored_list_of_format_keys.append("GT")
            set_of_format_keys.discard('GT')
        anchored_list_of_format_keys.extend(set_of_format_keys)
        return anchored_list_of_format_keys

    def merge_format_keys(self, other_vcf_line):
        """ Storing and merging of FORMAT keys of a VCF line.
        :param: other_vcf_line: the other VCF line to merge with
        """
        merged_format_keys = set()
        this_keys = self.format.split(":")
        other_keys = other_vcf_line.format.split(":")
        for this_key in this_keys:
            merged_format_keys.add(this_key)
        for other_key in other_keys:
            merged_format_keys.add(other_key)
        list_of_merged_format_key = self.order_format_keys(merged_format_keys)
        self.format = ":".join(list_of_merged_format_key)
        other_vcf_line.format = ":".join(list_of_merged_format_key)

    def combine_format_values_by_sample(self, format_tag_and_values_per_sample, list_of_sample_names):
        """ Creates a partial vcf data line of sample format values.
        :param format_tag_and_values_per_sample: nested dictionary {sample_name: {format_tag:formatvalue}}.
        :param list_of_sample_names: list of sample names
        :return: sample_format_values_string: formatted string (in the VCF file, this would be the tab-separated values under the sample name)
        """
        sample_format_value_tokens = []
        # Creates the list of FORMAT keys so we can get its corresponding value later
        set_of_format_keys = {key for sample in format_tag_and_values_per_sample for key in
                              format_tag_and_values_per_sample[sample]}
        list_of_format_key = self.order_format_keys(set_of_format_keys)
        # Generate string. For present samples, get its format value. For missing samples, populate with a missing value.
        for sample in list_of_sample_names:
            if sample in format_tag_and_values_per_sample:
                format_value_list = []
                for key in list_of_format_key:
                    format_value_list.append(format_tag_and_values_per_sample.get(sample, '.').get(key,
                                                                                                '.'))  # adds missing values if not found
                sample_format_value_tokens.append(":".join(format_value_list))
            else:
                sample_format_value_tokens.append(':'.join(['.' for key in list_of_format_key] or ['.']))
        self.format_values_by_sample_string = '\t'.join(sample_format_value_tokens)
        return self.format_values_by_sample_string
    # functions responsible for INFO are below
    def convert_info_list_to_dict(self):
        """ This converts list of INFO fields in a VCF line to a dictionary. This will be useful when merging fields of a VCF line.
        :return: info_dict: converted dictionary of INFO fileds
        """
        info_dict = {}
        for i in self.info:
            tokens = i.split(";")
            for token in tokens:
                info_key, info_value = token.split("=")
                info_dict[info_key] = info_value
        # NOTE: info_dict is a contains same info as vcf_value except the keys are converted
        return info_dict

    def merge_info_dicts(self, other_vcf_line):
        """ Merges and stores the INFO dictionaries for the INFO field of a VCF line.
        :param: other_vcf_line
        """
        merged_info_dict = {}
        for key in self.info_dict.keys() | other_vcf_line.info_dict.keys():
            self.value_info_dict = self.info_dict.get(key)
            other_vcf_line.value_info_dict = other_vcf_line.info_dict.get(key)
            if self.value_info_dict == other_vcf_line.value_info_dict:
                merged_info_dict[key] = self.value_info_dict
            else:
                if self.value_info_dict is None:
                    merged_info_dict[key] = other_vcf_line.value_info_dict
                elif other_vcf_line.value_info_dict is None:
                    merged_info_dict[key] = self.value_info_dict
                else:
                    merged_info_dict[key] = f"{self.value_info_dict},{other_vcf_line.value_info_dict}"
        # Store merged info dict for this VCF line and the other VCF line.
        self.info_dict = merged_info_dict
        other_vcf_line.info_dict = merged_info_dict

    def format_info_string(self):
        """ Creates a formatted INFO string using the INFO dictionary. Anchors ID to start of the string.
        :return: info_string: formatted INFO string for use in VCF line
        """
        info_parts = []
        if "ID" in self.info_dict:
            info_parts.append(f"ID={self.info_dict["ID"]}")
        for key, value in self.info_dict.items():
            if key != "ID":
                info_parts.append(f"{key}={value}")
        self.info_string = ";".join(info_parts)
        return self.info_string

    # MERGE OR KEEP below
    def merge(self, other_vcf_line, list_of_sample_names):
        """ Merging of the fields of a VCF line (ID, ALT, FILTER, INFO, FORMAT, FORMATvalues).
        :param: other_vcf_line : other VCF line to merge with
        :param: list_of_sample_names: list of sample names to help with creating format values by sample
        """
        # Merging INFO fields
        self.info_dict =  self.convert_info_list_to_dict()
        other_vcf_line.info_dict = other_vcf_line.convert_info_list_to_dict()
        self.merge_info_dicts(other_vcf_line)
        self.info_string  = self.format_info_string().rstrip(";")
        # merging FORMAT keys
        self.merge_format_keys(other_vcf_line)
        # merging FORMAT values
        merged_format_dict = self.format_dict | other_vcf_line.format_dict
        self.format_dict = merged_format_dict
        other_vcf_line.format_dict = merged_format_dict
        self.format_values_by_sample_string = self.combine_format_values_by_sample(self.format_dict, list_of_sample_names)
        other_vcf_line.format_values_by_sample_string = other_vcf_line.combine_format_values_by_sample(other_vcf_line.format_dict, list_of_sample_names)
        return (self.chrom,
                self.pos,
                self.merge_and_add(self.id, other_vcf_line.id, ";"),
                self.ref,
                self.merge_and_add(self.alt, other_vcf_line.alt, ","),
                self.qual,
                self.merge_and_add(self.filter, other_vcf_line.filter, ";"),
                self.info_string,
                self.format,
                self.format_values_by_sample_string
                )

    def keep(self, list_of_sample_names):
        self.format_values_by_sample_string = self.combine_format_values_by_sample(self.format_dict, list_of_sample_names)
        return (self.chrom,
                self.pos,
                self.id,
                self.ref,
                self.alt,
                self.qual,
                self.filter,
                self.info_string,
                self.format,
                self.format_values_by_sample_string
                )
