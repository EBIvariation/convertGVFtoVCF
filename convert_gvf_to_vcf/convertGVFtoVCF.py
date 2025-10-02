import argparse
import os

from Bio import SeqIO

# setting up paths to useful directories
convert_gvf_to_vcf_folder = os.path.dirname(__file__)
etc_folder = os.path.join(convert_gvf_to_vcf_folder, 'etc')

def read_file(prefix, header_type):
    """Reads in {reserved/sv}{INFO/FORMAT}keys.tsv files and returns the dictionary where the key is the KEYID
    (usually column 1) and the value is a list of file tokens
    :param prefix: prefix of files to read i.e. sv or reserved
    :param header_type: type of header file to read i.e. INFO or FORMAT
    :return: file lines: a dictionary where the key is the KEYID and the value is a list of file tokens
    """
    file_lines = {}
    keys_tsv_file = os.path.join(etc_folder, f'{prefix}{header_type}keys.tsv')
    try:
        with open(keys_tsv_file) as keys_file:
            next(keys_file) # Skip the header
            for line in keys_file:
                file_tokens = line.rstrip().split("\t")
                key_id = file_tokens[0]
                number_of_tokens = len(file_tokens)
                if number_of_tokens >= 2:
                    for token in range(number_of_tokens):
                        value_to_add = file_tokens[token]
                        file_lines.setdefault(key_id, []).append(value_to_add)
    except FileNotFoundError:
        print(f'File not found: {keys_tsv_file}')
    return file_lines

def generate_custom_structured_metainformation_line(vcf_key, vcf_key_id, vcf_key_number, vcf_key_type, vcf_key_description,
                                                    optional_extra_fields=None):
    """ Generates a custom structured meta-information line for INFO/FILTER/FORMAT/ALT
    :param vcf_key: required field INFO, FILTER, FORMAT, ALT
    :param vcf_key_id: required field for structured lines ID
    :param vcf_key_number: The number of values that can be included or special character: A or R or G or .
    :param vcf_key_type: Values are Integer, Float, Character, String
    :param vcf_key_description: Description
    :param optional_extra_fields: an optional field, dictionary of custom fields and their values
    :return: custom_structured_string
    """
    extra_keys_kv_lines = []
    if optional_extra_fields:
        for extra_field in optional_extra_fields:
            kv_line = "," + extra_field + "=" + '"' + optional_extra_fields[extra_field] + '"'
            extra_keys_kv_lines.append(kv_line)
    vcf_key_extra_keys = ''.join(extra_keys_kv_lines)
    custom_structured_string = f'##{vcf_key}=<ID="{vcf_key_id}",Number="{vcf_key_number}",Type="{vcf_key_type}",Description="{vcf_key_description}"{vcf_key_extra_keys}>'
    return custom_structured_string

def generate_all_possible_standard_structured_lines(header_type):
    """ Generates a fictionary of all possible standard structured lines for INFO/FILTER/FORMAT/ALT
    :param header_type: type of header file to read i.e. ALT, FILTER, INFO or FORMAT
    :return: dictionary of all possible standard structured lines keys for the header type
    """
    all_possible_lines = {}
    prefix = {}
    prefix["reserved"] = True
    prefix["sv"] = True
    if header_type == 'ALT':
        prefix["reserved"] = False
    if prefix["reserved"]:
        reserved_key = read_file("reserved", header_type)
        for r_key in reserved_key:
            key_id = reserved_key[r_key][0]
            number = reserved_key[r_key][1]
            type_for_key = reserved_key[r_key][2]
            description = reserved_key[r_key][3]
            reserved_string = f'##{header_type}=<ID={key_id},Number={number},Type={type_for_key},Description="{description}">'
            all_possible_lines[key_id] = reserved_string
    if prefix['sv']:
        sv_key = read_file("sv", header_type)
        for s_key in sv_key:
            sv_key_id = sv_key[s_key][0]
            sv_line = sv_key[s_key][1]
            all_possible_lines[sv_key_id] = sv_line
    return all_possible_lines


def generate_standard_structured_metainformation_line(vcf_key_id, standard_lines_for_vcf_key, all_possible_lines):
    """Generates a list of standard structured metainformation lines.
    :param vcf_key_id: VCF tag key id
    :param standard_lines_for_vcf_key: lines_standard_NAME i.e a list of standard lines for this VCF file regarding INFO or ALT or FILTER or FORMAT
    :param all_possible_lines: all_possible_NAME_lines i.e. list of all possible lines for INFO or ALT or FILTER or FORMAT
    :return: standard_lines_for_vcf_key: a dictionary
    """
    standard_structured_line = all_possible_lines[vcf_key_id]
    standard_lines_for_vcf_key.append(standard_structured_line)
    return standard_lines_for_vcf_key


def generate_custom_unstructured_metainformation_line(vcf_unstructured_key, vcf_unstructured_value, lines_custom_unstructured):
    """ Generates a formatted unstructured metainformation line using a custom key value pair. This is stored in the list called lines_custom_unstructured.
    :param lines_custom_unstructured: list to store custom unstructured metainformation lines
    :param vcf_unstructured_key: key for custom unstructured metainformation line
    :param vcf_unstructured_value: value for custom unstructured metainformation line
    :return: custom_unstructured_string
    """
    custom_unstructured_string = f"##{vcf_unstructured_key}={vcf_unstructured_value}"
    lines_custom_unstructured.append(custom_unstructured_string)
    return custom_unstructured_string

def read_info_attributes(info_attributes_file):
    """ Read in the file containing specific INFO attributes.
    :param info_attributes_file: A file containing the specific attributes
    :return: attribute_dict: A dictionary of key id and the list of attribute tokens
    """
    attribute_dict = {}  # dictionary of dgva specific INFO attributes
    with open(info_attributes_file) as open_file:
        next(open_file)
        for line in open_file:
            attribute_tokens = line.rstrip().split("\t")
            key = attribute_tokens[0]
            attribute_dict[key] = attribute_tokens
    return attribute_dict

def read_sequence_ontology_symbolic_allele(so_symbolic_allele_file):
    """ Read in the file containing sequence ontology symbolic allele and returns a dictionary.
    :param: so_symbolic_allele_file - the file of sequence ontology symbolic alleles.
    :return: symbolic allele dictionary - symbolic alleles as key and list of variant types as the value.
    """
    symbolic_allele_dict = {}
    with open(so_symbolic_allele_file) as so_symbolic_allele:
        next(so_symbolic_allele)
        for line in so_symbolic_allele:
            allele_tokens = line.rstrip().split("\t")
            symb_allele = allele_tokens[0]
            sequence_ontology_id = allele_tokens[2]
            name = allele_tokens[3]
            description = allele_tokens[4]
            symbolic_allele_dict.setdefault(name, []).append(sequence_ontology_id)
            symbolic_allele_dict.setdefault(name, []).append(symb_allele)
            symbolic_allele_dict.setdefault(name, []).append(description)
    return symbolic_allele_dict

def extract_reference_allele(fasta_file, chromosome_name, position, end):
    """ Extracts the reference allele from the assembly.
    :param fasta_file: FASTA file of the assembly
    :param chromosome_name: name of the sequence
    :param position: position
    :param end: end position
    :return: reference_allele: base found at this chromosome_name at this position within this fasta_file
    """
    # using .index instead of .todict for memory efficiency:  https://biopython.org/docs/1.76/api/Bio.SeqIO.html#input-multiple-records
    records_dictionary = SeqIO.index(fasta_file, "fasta")
    zero_indexed_position = position - 1 # minus one because zero indexed
    zero_indexed_end = end - 1
    reference_allele = ""
    for position in range(zero_indexed_position, zero_indexed_end):
        reference_allele = reference_allele + records_dictionary[chromosome_name].seq[position]
    records_dictionary.close()
    return reference_allele

def get_gvf_attributes(column9_of_gvf):
    """Get a dictionary of GVF attributes
    :param column9_of_gvf:  column - the final column of the GVF file
    :return: gvf_attribute_dictionary: a dictionary of attribute keys and their values
    """
    gvf_attribute_dictionary = {} # attribute key => value
    # parse by semicolon this creates attribute
    # parse by equals sign this creates tag-values, if the value is a comma, create a list
    attributes_in_gvf_line = column9_of_gvf.split(";")
    for attribute in attributes_in_gvf_line:
        attribute_key, attribute_value = attribute.split("=")
        if "," in attribute_value:
            attribute_value_list = attribute_value.split(",")
            gvf_attribute_dictionary[attribute_key] = attribute_value_list
        else:
            gvf_attribute_dictionary[attribute_key] = attribute_value
    return gvf_attribute_dictionary


# CAVEATS: 1) assume sample_name is present in the GVF file. If absent consider adding UnknownSample1, UnknownSample2 etc.
def convert_gvf_attributes_to_vcf_values(column9_of_gvf,
                                         dgva_attribute_dict,
                                         gvf_attribute_dict,
                                         lines_custom_structured,
                                         standard_lines_dictionary,
                                         all_possible_lines_dictionary):
    gvf_attribute_dictionary = get_gvf_attributes(column9_of_gvf)
    vcf_vals = {}
    catching_for_review = []
    # created a rough guide to attributes_for_custom_structured_metainformation in dgvaINFOattributes.tsv = this probably should be refined at a later date
    # TODO: edit dgvaINFOattributes.tsv i.e. replace unknown placeholders '.' with the actual answer, provide a more informative description
    for attrib_key in gvf_attribute_dictionary:
        # if dgva specific key, create custom string otherwise do standard
        if attrib_key in dgva_attribute_dict:
            lines_custom_structured.append(
                generate_custom_structured_metainformation_line(vcf_key="INFO", vcf_key_id=attrib_key,
                                                                vcf_key_number=dgva_attribute_dict[attrib_key][1],
                                                                vcf_key_type=dgva_attribute_dict[attrib_key][2],
                                                                vcf_key_description=dgva_attribute_dict[attrib_key][3],
                                                                optional_extra_fields=None)
            )
            vcf_vals[attrib_key]=gvf_attribute_dictionary[attrib_key]
        elif attrib_key == "allele_count":
            #generate_standard_structured_metainformation_line("INFO", "AC", lines_standard_ALT, lines_standard_INFO, lines_standard_FILTER, lines_standard_FORMAT, all_possible_ALT_lines, all_possible_INFO_lines, all_possible_FILTER_lines, all_possible_FORMAT_lines)
            lines_standard_info = generate_standard_structured_metainformation_line("AC", standard_lines_dictionary["INFO"], all_possible_lines_dictionary["INFO"])
        elif attrib_key == "allele_frequency":
            lines_standard_info = generate_standard_structured_metainformation_line("AF", standard_lines_dictionary["INFO"], all_possible_lines_dictionary["INFO"])
        elif attrib_key == "ciend":
            lines_standard_info = generate_standard_structured_metainformation_line("CIEND", standard_lines_dictionary["INFO"], all_possible_lines_dictionary["INFO"])
        elif attrib_key == "copy_number":
            lines_standard_info = generate_standard_structured_metainformation_line("CN", standard_lines_dictionary["INFO"], all_possible_lines_dictionary["INFO"])
        elif attrib_key == "insertion_length":
            lines_standard_info = generate_standard_structured_metainformation_line("SVLEN", standard_lines_dictionary["INFO"],
                                                                                    all_possible_lines_dictionary["INFO"])
        elif attrib_key == "mate_id":
            lines_standard_info = generate_standard_structured_metainformation_line("MATEID", standard_lines_dictionary["INFO"],
                                                                                    all_possible_lines_dictionary["INFO"])
        elif attrib_key == "sample_name":
            #sample_names.append(sample_names)
            pass
        # GVF keys (not dgva specific)
        elif attrib_key == "ID":
            pass
        elif attrib_key == "Variant_seq":
            pass
        elif attrib_key == "Reference_seq":
            pass
        elif (attrib_key == "Alias" or attrib_key == "Variant_effect" or attrib_key == "Variant_codon" or
              attrib_key == "Reference_codon" or attrib_key == "Variant_aa" or attrib_key == "Reference_aa" or
              attrib_key == "breakpoint_detail" or attrib_key == "Sequence_context"):
            lines_custom_structured.append(
                generate_custom_structured_metainformation_line(
                    vcf_key="INFO", vcf_key_id=attrib_key,
                    vcf_key_number=gvf_attribute_dict[attrib_key][1],
                    vcf_key_type=gvf_attribute_dict[attrib_key][2],
                    vcf_key_description=gvf_attribute_dict[attrib_key][3],
                    optional_extra_fields=None)
            )
        elif attrib_key == "Dbxref":
            # custom info tag + pase and add to id?
            pass
        elif attrib_key == "Variant_reads":
            # reserved info/format key, AD/AC
            pass
        elif attrib_key == "Total_reads":
            # reserved info key, DP
            pass
        elif attrib_key == "Variant_freq":
            # reserve info tag, AF
            pass
        elif attrib_key == "Zygosity":
            # format and GT tag
            pass
        elif attrib_key == "Genotype":
            # GT
            pass
        elif attrib_key == "Phased":
            # GT or FORMAT PS
            pass
        elif attrib_key == "Start_range":
            # either custom info tag or CIPOS or CIEND, may need imprecise
            pass
        elif attrib_key == "End_range":
            # either custom info tag or CIEND, may need imprecise
            pass
        elif attrib_key == "Breakpoint_range":
            # either custom info tag or CIPOS, CIEND, may need imprecise
            pass
        elif attrib_key == "Individual":
            # sampl name for each column
            pass
        else:
            print("catching these attribute keys for review at a later date", attrib_key)
            catching_for_review.append(attrib_key)
    #print("dictionary", gvf_attribute_dictionary)
    # print("vcf_vals", vcf_vals)
    return gvf_attribute_dictionary

# step 7c
class GvfFeatureline:
    def __init__(self, seqid, source, so_type, start, end, score, strand, phase, attributes):
        """ Initialise GvfFeature line
        :param seqid: sequence ID i.e. chromosome number
        :param source: source i.e. DGVa
        :param so_type: sequence ontology structural variant type
        :param start: start position
        :param end: end position
        :param score: score
        :param strand: strandedness
        :param phase: phase
        :param attributes: attribute key values
        """
        self.seqid = seqid
        self.source = source
        self.feature_type = so_type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = attributes

    def __str__(self):
        """
        Helper to print variables of the GVF feature line
        :return:line_to_print
        """
        line_to_print = self.seqid + "\t" + self.source + "\t" + self.feature_type + "\t" + self.start + "\t" + self.end + "\t" + self.score +"\t" + self.strand + "\t" + self.phase + "\t" + self.attributes
        return line_to_print

# step 7
def read_in_gvf_file(gvf_input):
    """ Reads in the user provided GVF file.
    :param gvf_input: arguments.gvf_input
    :return: gvf_pragmas, gvf_non_essential, gvf_lines_obj_list
    """
    gvf_pragmas = []  # list of pragma lines starting with: ##
    gvf_non_essential = []  # list of non-essential lines starting with: #
    features = []
    gvf_lines_obj_list = []  # list of objects when reading in gvf files, one object represents a gvf line
    print("Reading in the following GVF input: " + gvf_input)
    with open(gvf_input) as gvf_file:
        for line in gvf_file:
            if line.startswith("##"):
                gvf_pragmas.append(line.rstrip())
            elif line.startswith("#"):
                gvf_non_essential.append(line.rstrip())
            else:
                features.append(line.rstrip())
    for feature in features:
        f_list = feature.split("\t")
        line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6], f_list[7], f_list[8])
        gvf_lines_obj_list.append(line_object)
    return gvf_pragmas, gvf_non_essential, gvf_lines_obj_list

#step 8
#TODO: ID this can be a semi-colon separated list or a '.' (if no value = '.'; one value = value; more than one = value;value)
class VcfLine:
    def __init__(self, gvf_feature_line_object,
                 dgva_attribute_dict,
                 gvf_attribute_dict,
                 symbolic_allele_dictionary,
                 assembly_file,
                 lines_custom_structured,
                 standard_lines_dictionary,
                 all_possible_lines_dictionary):

        # ATTRIBUTES
        self.vcf_value = convert_gvf_attributes_to_vcf_values(gvf_feature_line_object.attributes,
                                                              dgva_attribute_dict,
                                                              gvf_attribute_dict,
                                                              lines_custom_structured,
                                                              standard_lines_dictionary,
                                                              all_possible_lines_dictionary)

        self.assembly = assembly_file
        self.symbolic_allele_dictionary = symbolic_allele_dictionary

        # GVF
        self.source = gvf_feature_line_object.source
        self.so_type = gvf_feature_line_object.feature_type #currently column 3 of gvf, but could be an attribute so perhapsVCF: INFO or FORMAT?
        self.end = int(gvf_feature_line_object.end)
        self.phase = gvf_feature_line_object.phase # this is always a placeholder'.'

        # VCF DATALINE
        self.chrom = gvf_feature_line_object.seqid
        self.pos = int(gvf_feature_line_object.start)
        self.id = self.vcf_value["ID"]  # attributes: ID
        self.length = self.end - self.pos
        self.qual = gvf_feature_line_object.score # see EVA-3879: this is always '.'
        self.filter = "." # this is always a placeholder '.'; perhaps could add s50.

        # INFO
        #TODO: specific SV info keys populated from gvf_feature_line
        self.key = self.chrom + "_" + str(self.pos)
        self.info = [] # TODO: add info field for self.info
        # calculated last
        self.ref = self.get_ref()
        self.alt = self.get_alt(standard_lines_dictionary, all_possible_lines_dictionary)

        self.sample_name = self.vcf_value["sample_name"] # this should be each samples names format value # sample names needs to be populated in attributes
        # # higher priority
        self.format = "pending" #TODO: set this in convertgvfattributes


        # # each item in the list exclude_from_info has its own place in the VCF file, so not part of info
        # exclude_from_info = ["ID", # done above
        #                      "Variant_seq", # done above
        #                      "Reference_Seq", # done above
        #                      "Genotype",
        #                      "Phased",
        #                      "sample_name",
        #                      "Zygosity"]
        # for k in self.attributes.keys():
        #     if k not in exclude_from_info:
        #         self.info = self.info + k + "=" + str(self.attributes[k]) + ";"
        #     elif k =="ID":
        #         pass
        #     elif k == "Variant_seq":
        #         #TODO: ADD TO ALT
        #         pass
        #     elif k == "Reference_Seq":
        #         # TODO: ADD TO REF
        #         pass
        #     elif k == "Genotype":
        #         # TODO: ADD TO GT (FORMAT)
        #         pass
        #     elif k == "Phased":
        #         # TODO: ADD TO GT
        #         pass
        #     elif k == "sample_name":
        #         # ADD TO LIST: sample_name
        #         pass
        #     elif k == "Zygosity":
        #         # add to format, gt flag
        #         pass
        #     else:
        #         pass

    def add_padded_base(self, ref, alt, placed_before : bool):
        """ Adds padded base to REF and ALT allele
        :param ref: reference allele
        :param alt: alt allele
        :param placed_before: padded base is placed before ref or alt True or False
        :return: (padded_base, self.pos, self.ref, self.alt)
        """
        if placed_before:
            padded_base_pos = self.pos - 1
            self.pos = padded_base_pos
            padded_base = extract_reference_allele(self.assembly, self.chrom, self.pos, self.end)
            ref = padded_base + ref
            if alt == ".":
                alt = padded_base
            else:
                alt = padded_base + alt
        elif not placed_before:
            padded_base_pos = self.pos + 1
            new_end = self.end + 1
            padded_base = extract_reference_allele(self.assembly, self.chrom, padded_base_pos, new_end)
            ref = ref + padded_base
            if alt == ".":
                alt = padded_base
            else:
                alt = alt + padded_base
        else:
            print("WARNING: Variable placed_before unknown: " + str(placed_before))
            padded_base = None
        return padded_base, self.pos, ref, alt

    def build_iupac_ambiguity_code(self):
        """ Builds dictionary for the iupac ambiguity code
        :return: iupac_ambiguity_dictionary: iupac code as key, list of values as value
        """
        # see PMID: 20202974 (Table 1) for the official list
        iupac_codes = ["R", "Y", "M", "K", "S", "D", "W", "H", "B", "V", "D", "N"]
        R = ["A", "G"]
        Y = ["C", "T"]
        M = ["A", "C"]
        K = ["G", "T"]
        S = ["C", "G"]
        W = ["A", "T"]
        H = ["A", "C", "T"]
        B = ["C", "G", "T"]
        V = ["A", "C", "G"]
        D = ["A", "G", "T"]
        N = ["A", "C", "G", "T"]
        iupac_values = [R, Y, M, K, S, D, W, H, B, V, D, N]
        iupac_ambiguity_dictionary = dict(zip(iupac_codes, iupac_values))
        return iupac_ambiguity_dictionary

    def convert_iupac_ambiguity_code(self, iupac_ambiguity_dictionary, ref_to_convert):
        """ Converts the REF allele if it contains IUPAC ambiguity cod
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

    def check_ref(self, ref_allele_to_be_checked):
        """ Checks whether a reference allele meets the requirements of the VCF specification
        :param ref_allele_to_be_checked: reference allele to check
        :return: checked_reference_allele: reference allele that meets the requirements of the VCF specification"""
        if isinstance(ref_allele_to_be_checked, str):
            if not all(bases in ref_allele_to_be_checked for bases in ["A", "C", "G", "T", "N"]):
                iupac_ambiguity_dictionary = self.build_iupac_ambiguity_code()
                checked_reference_allele = self.convert_iupac_ambiguity_code(iupac_ambiguity_dictionary, ref_allele_to_be_checked)
            else:
                checked_reference_allele = ref_allele_to_be_checked
        else:
            print("WARNING: Ref allele must be a string. Setting to missing value.", ref_allele_to_be_checked)
            checked_reference_allele = "."
        return checked_reference_allele

    def get_ref(self):
        """ Gets the reference allele from attributes column or if not found, returns "."
        :return: reference allele
        """
        if "Reference_seq" in self.vcf_value.keys():
            reference_allele = self.vcf_value["Reference_seq"]
        else:
            if self.assembly:
                reference_allele = extract_reference_allele(self.assembly, self.chrom, self.pos, self.end)
            else:
                print("WARNING: No reference provided. Placeholder inserted for Reference allele.")
                reference_allele = "."
        if reference_allele != ".":
            reference_allele = self.check_ref(reference_allele)
        return reference_allele


    def generate_symbolic_allele(self, standard_lines_dictionary, all_possible_lines_dictionary):
        """ Generates the symbolic allele and stores the corresponding metainformation lines. Also determines if variant is precise or imprecise.
        :param lines_standard_alt: stores ALT metainformation lines
        :param lines_standard_info: stores INFO metainformation lines
        :param all_possible_alt_lines: list of all possible ALT lines
        :param all_possible_info_lines: list of all possible INFO lines
        :return: symbolic_allele, self.info, lines_standard_ALT, lines_standard_INFO
        """
        symbolic_allele_id = self.symbolic_allele_dictionary[self.so_type][1]
        symbolic_allele = f'<{symbolic_allele_id}>'

        lines_standard_alt = standard_lines_dictionary["ALT"]
        lines_standard_info = standard_lines_dictionary["INFO"]
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
                "Cannot identify symbolic allele"
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

    def get_alt(self, standard_lines_dictionary, all_possible_lines_dictionary):
        """ Gets the ALT allele for the VCF file
        :param lines_standard_alt: store ALT lines
        :param lines_standard_info: store NFO lines
        :param all_possible_alt_lines: dictionary of all possible ALT lines
        :param all_possible_info_lines: dictionary of all possible INFO lines
        :return: symbolic_allele, self.info, lines_standard_ALT, lines_standard_INFO
        """
        lines_standard_alt = standard_lines_dictionary["ALT"]
        lines_standard_info = standard_lines_dictionary["INFO"]
        all_possible_alt_lines = all_possible_lines_dictionary["ALT"]
        all_possible_info_lines = all_possible_lines_dictionary["INFO"]

        if any(base in self.vcf_value["Variant_seq"] for base in ["A", "C", "G", "T", "N"]):
            alterative_allele = self.vcf_value["Variant_seq"]
        elif self.vcf_value["Variant_seq"] == '.':
            symbolic_allele, self.info, lines_standard_alt, lines_standard_info = self.generate_symbolic_allele(standard_lines_dictionary, all_possible_lines_dictionary)
            if symbolic_allele is None:
                alterative_allele = "."
            elif (self.vcf_value["Variant_seq"] == "." or self.vcf_value["Variant_seq"] == "-") and symbolic_allele is not None:
                alterative_allele = symbolic_allele
                # add padded bases
                if self.pos == 1:
                    #print("pos, ref, alt",self.pos,self.ref, alterative_allele)
                    padded_base, self.pos, self.ref, self.alt = self.add_padded_base(self.ref, alterative_allele, False)
                    self.ref = self.check_ref(self.ref)
                else:
                    #print("pos, ref, alt", self.pos,self.ref, alterative_allele)
                    padded_base, self.pos, self.ref, self.alt = self.add_padded_base(self.ref, alterative_allele, True)
                    self.ref = self.check_ref(self.ref)
            else:
                alterative_allele = "."
                print("Cannot identify symbolic allele. Variant type is not supported.")
        else:
            alterative_allele = "."
            print("Could not determine the alterative allele.")
        return alterative_allele

    def __str__(self):
        string_to_return = '\t'.join((self.chrom, self.pos, self.key, self.qual, self.filter, self.info, self.source, self.phase, self.end, self.so_type, self.sample_name, self.format))
        return string_to_return

#step 9 using custom unstructured meta-information line = generate_custom_unstructured_metainfomation_line
def generate_vcf_metainformation(lines_custom_unstructured, gvf_pragmas, gvf_non_essential, list_of_vcf_objects,
                                 standard_lines_dictionary):
    """ Generates a list of metainformation lines for the VCF header
    :param lines_custom_unstructured: a list of formatted unstructured metainformation lines using a custom key value pair
    :param gvf_pragmas: list of gvf pragmas to convert
    :param gvf_non_essential: list of non-essential gvf pragmas to convert
    :param list_of_vcf_objects: list of vcf objects
    :param lines_standard_alt: list of ALT lines
    :param lines_standard_info: list of INFO lines
    :param lines_standard_filter: list of FILTER lines
    :param lines_standard_format: list of FORMAT lines
    :return: unique_pragmas_to_add, sample_names: a list of pragmas (this list contains no duplicates), list of sample names
    """
    pragmas_to_add = []
    unique_pragmas_to_add = []
    sample_names = []

    unique_alt_lines_to_add = []
    unique_info_lines_to_add = []
    unique_filter_lines_to_add = []
    unique_format_lines_to_add = []
    # MANDATORY: file format for VCF
    pragma_fileformat = generate_custom_unstructured_metainformation_line("fileformat", "VCFv4.4", lines_custom_unstructured)
    pragmas_to_add.append(pragma_fileformat)
    #Go through essential pragmas
    #TODO: list of pragmas to add:reference=file, contig, phasing,INFO#
    for pragma in gvf_pragmas:
        # file date
        if pragma.startswith("##file-date"):
            date = pragma.split(" ")[1].replace("-", "")
            pragma_filedate = generate_custom_unstructured_metainformation_line("fileDate", date, lines_custom_unstructured)
            pragmas_to_add.append(pragma_filedate)
        # source
        for vcf_obj in list_of_vcf_objects:
            pragma_source = generate_custom_unstructured_metainformation_line("source", vcf_obj.source, lines_custom_unstructured)
            pragmas_to_add.append(pragma_source)
        # reference #TODO: add this
        # contig: recommended, see section 1.4.7 of VCF specification # TODO: add this
        # phasing # not required
        if pragma.startswith("##gff-version"):
            gff_version_number = pragma.split(" ")[1]
            pragma_gff_version_number = generate_custom_unstructured_metainformation_line("gff-version", gff_version_number, lines_custom_unstructured)
            pragmas_to_add.append(pragma_gff_version_number)
        elif pragma.startswith("##gvf-version"):
            gvf_version_number = pragma.split(" ")[1]
            pragma_gvf_version_number = generate_custom_unstructured_metainformation_line("gvf-version", gvf_version_number, lines_custom_unstructured)
            pragmas_to_add.append(pragma_gvf_version_number)
        elif pragma.startswith("##species"):
            species_value = pragma.split(" ")[1]
            pragma_species_value = generate_custom_unstructured_metainformation_line("species", species_value, lines_custom_unstructured)
            pragmas_to_add.append(pragma_species_value)
        elif pragma.startswith("##genome-build"):
            genome_build = pragma.split("genome-build ")[1]
            pragma_genome_build = generate_custom_unstructured_metainformation_line("genome-build", genome_build, lines_custom_unstructured)
            pragmas_to_add.append(pragma_genome_build)
        else:
            pass
    # Go through non-essential pragmas
    for non_essential_pragma in gvf_non_essential:
        if non_essential_pragma.startswith("#Study_accession"):
            study_accession = non_essential_pragma.split(": ")[1]
            non_essential_pragma_study_accession = generate_custom_unstructured_metainformation_line("Study_accession", study_accession, lines_custom_unstructured)
            pragmas_to_add.append(non_essential_pragma_study_accession)
        elif non_essential_pragma.startswith("#Study_type"):
            study_type = non_essential_pragma.split(": ")[1]
            non_essential_pragma_study_type = generate_custom_unstructured_metainformation_line("Study_type", study_type, lines_custom_unstructured)
            pragmas_to_add.append(non_essential_pragma_study_type)
        elif non_essential_pragma.startswith("#Display_name"):
            display_name = non_essential_pragma.split(": ")[1]
            non_essential_pragma_display_name = generate_custom_unstructured_metainformation_line("Display_name", display_name, lines_custom_unstructured)
            pragmas_to_add.append(non_essential_pragma_display_name)
        elif non_essential_pragma.startswith("#Publication"):
            publication = non_essential_pragma.split(": ")[1]
            non_essential_pragma_publication = generate_custom_unstructured_metainformation_line("Publication", publication, lines_custom_unstructured)
            pragmas_to_add.append(non_essential_pragma_publication)
        elif non_essential_pragma.startswith("#Study"):
            study = non_essential_pragma.split(": ")[1]
            non_essential_pragma_study = generate_custom_unstructured_metainformation_line("Study", study, lines_custom_unstructured)
            pragmas_to_add.append(non_essential_pragma_study)
        elif non_essential_pragma.startswith("#Assembly_name"):
            assembly_name = non_essential_pragma.split(": ")[1]
            non_essential_pragma_assembly_name = generate_custom_unstructured_metainformation_line("Assembly_name", assembly_name, lines_custom_unstructured)
            pragmas_to_add.append(non_essential_pragma_assembly_name)
        elif non_essential_pragma.startswith("#subject"):
            subject = non_essential_pragma.split(": ")[1]
            non_essential_pragma_subject = generate_custom_unstructured_metainformation_line("subject", subject, lines_custom_unstructured)
            pragmas_to_add.append(non_essential_pragma_subject)
        elif non_essential_pragma.startswith("#sample"):
            sample_information = non_essential_pragma.split(": ")[1]
            non_essential_pragma_sample = generate_custom_unstructured_metainformation_line("sample", sample_information,
                                                                                            lines_custom_unstructured)
            pragmas_to_add.append(non_essential_pragma_sample)
            list_of_sample_information = sample_information.split(";")
            for sample_info in list_of_sample_information:
                if sample_info.startswith("sample_name"):
                    sample_name = sample_info.split("=")[1]
                    sample_names.append(sample_name)
        else:
            print("Skipping unknown non-essential GVF pragma:", non_essential_pragma)

    print("Total number of samples in this VCF: ", len(sample_names))

    for pragma in pragmas_to_add:
        if pragma not in unique_pragmas_to_add:
            unique_pragmas_to_add.append(pragma)
    for alt_line in standard_lines_dictionary["ALT"]:
        if alt_line not in unique_alt_lines_to_add:
            unique_alt_lines_to_add.append(alt_line)
    for info_line in standard_lines_dictionary["INFO"]:
        if info_line not in unique_info_lines_to_add:
            unique_info_lines_to_add.append(info_line)
    for filter_line in standard_lines_dictionary["FILTER"]:
        if filter_line not in unique_filter_lines_to_add:
            unique_filter_lines_to_add.append(filter_line)
    for format_line in standard_lines_dictionary["FORMAT"]:
        if format_line not in unique_format_lines_to_add:
            unique_format_lines_to_add.append(format_line)
    return unique_pragmas_to_add, sample_names, unique_alt_lines_to_add, unique_info_lines_to_add, unique_filter_lines_to_add, unique_format_lines_to_add

# step 10
def generate_vcf_header_line(samples):
    """ Generates the VCF header line
    :param samples: list of samples, these will appear in the header line
    :return: vcf_header: a string
    """
    vcf_header_fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    for sample in samples:
        vcf_header_fields.append(sample)
    vcf_header = '\t'.join(vcf_header_fields)
    return vcf_header

def gvf_features_to_vcf_objects(gvf_lines_obj_list,
                                dgva_attribute_dict,
                                gvf_attribute_dict,
                                symbolic_allele_dictionary,
                                assembly_file,
                                lines_custom_structured,
                                standard_lines_dictionary,
                                all_possible_lines_dictionary):

    """ Creates VCF objects from GVF feature lines and stores the VCF objects.
    :param gvf_lines_obj_list: list of GVF feature line objects
    :param gvf_attribute_dict: dictionary of GVF INFO attributes
    :param symbolic_allele_dictionary: dictionary of symbolic alleles
    :param assembly_file: FASTA file to assembly
    :param dgva_attribute_dict: dictionary af DGVa specific INFO attributes
    :param lines_custom_structured: list to store custom structured metainformation lines
    :param lines_standard_alt: ALT lines for this VCF file
    :param lines_standard_info: INFO lines for this VCF file
    :param lines_standard_filter: FILTER lines for this VCF file
    :param lines_standard_format: FORMAT lines for this VCF file
    :param all_possible_lines_dictionary: dict of dictionaries. all possible ALT/INFO/FILTER/FORMAT lines
    :return: vcf_data_lines, list_of_vcf_objects: dictionary of lists and a list of VCF objects
    """
    vcf_data_lines = {}  # DICTIONARY OF LISTS
    list_of_vcf_objects = []

    # create a vcf object for every feature line in the GVF (1:1)
    # add the newly created vcf object to the vcf data line it belongs to
    # (1:many; key=chrom_pos; 1 key: many vcf objects)
    for gvf_featureline in gvf_lines_obj_list:
        vcf_object = VcfLine(gvf_featureline,
                             dgva_attribute_dict,
                             gvf_attribute_dict,
                             symbolic_allele_dictionary,
                             assembly_file,
                             lines_custom_structured,
                             standard_lines_dictionary,
                             all_possible_lines_dictionary)

        list_of_vcf_objects.append(vcf_object)
        if vcf_object.key in vcf_data_lines:
            vcf_data_lines[vcf_object.key].append(vcf_object)
        else:
            vcf_data_line_objects_list = []
            vcf_data_line_objects_list.append(vcf_object)
            vcf_data_lines[vcf_object.key] = vcf_data_line_objects_list
        # check the number of objects to see if they are merged
        # for key in vcf_data_lines.keys():
        #     vcf_obj_list = vcf_data_lines[key]
            # print("for ", key, " the number of vcf objects is: ", len(vcf_obj_list))
    return vcf_data_lines, list_of_vcf_objects


def populate_sample_formats(list_of_sample_names):
    """ Populates a dictionary using a list of sample names. Dictionary key is sample name, value is the sample's format value.
    :param list_of_sample_names: list of sample names
    :return:sample_name_format_value: dictionary of sample names => sample format value
    """
    sample_name_format_value = {}
    for sample in list_of_sample_names:
        sample_name_format_value[sample] = "sampleFORMAThere" #TODO: fill this in
    return sample_name_format_value

def format_sample_values(sample_name_format_value):
    """ Creates a partial vcf data line of sample format values.
    :param sample_name_format_value: dictionary of sample names => sample format value
    :return: sample_format_values_string: formatted string
    """
    sample_format_values_string = ""
    for key in sample_name_format_value:
        sample_format_values_string = sample_format_values_string + sample_name_format_value[key] + "\t"
    return sample_format_values_string

def format_vcf_datalines(list_of_vcf_objects, list_of_sample_names):
    """ Iterates through a list of VCF objects and sample names and formats them as a VCF dataline.
    :param list_of_vcf_objects: list of vcf objects
    :param list_of_sample_names: list of sample names
    :return: formatted_vcf_datalines: list of formatted vcf datalines
    """
    sample_name_format_value = populate_sample_formats(list_of_sample_names)
    sample_format_values_string = format_sample_values(sample_name_format_value)

    formatted_vcf_datalines = []
    for vcf_obj in list_of_vcf_objects:
        vcf_info_string = ";".join(vcf_obj.info)
        vcf_line = (f"{vcf_obj.chrom}\t"
                        f"{vcf_obj.pos}\t"
                        f"{vcf_obj.id}\t"
                        f"{vcf_obj.ref}\t" #TODO: should this always be empty
                        f"{vcf_obj.alt}\t" #TODO: should this always be empty
                        f"{vcf_obj.qual}\t" #TODO: should this always be empty
                        f"{vcf_obj.filter}\t" #TODO: should this always be empty
                        #f"{vcf_obj.info}\t"
                        f"{vcf_info_string}\t"
                        f"{vcf_obj.format}\t"
                        f"{sample_format_values_string}"
                        )
        formatted_vcf_datalines.append(vcf_line)
    return formatted_vcf_datalines

def main():
    print("Running the GVF to VCF converter")
    # step 1
    parser = argparse.ArgumentParser()
    parser.add_argument("gvf_input", help="GVF input file.")
    parser.add_argument("vcf_output", help="VCF output file.")
    parser.add_argument("-a", "--assembly", help="FASTA assembly file")
    args = parser.parse_args()
    print("The provided input file is: ", args.gvf_input)
    print("The provided output file is: ", args.vcf_output)
    if args.assembly:
        print("The provided assembly file is: ", args.assembly)

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

    # custom meta-information lines for this VCF file
    lines_custom_structured = []
    lines_custom_unstructured = []
    # standard structured meta-information lines for this VCF file
    lines_standard_alt = []
    lines_standard_info = []
    lines_standard_filter = []
    lines_standard_format = []
    # merging
    standard_lines_dictionary ={
        "ALT": lines_standard_alt,
        "INFO": lines_standard_info,
        "FILTER": lines_standard_filter,
        "FORMAT": lines_standard_format,
    }

    gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(args.gvf_input)
    dgva_info_attributes_file = os.path.join(etc_folder, 'dgvaINFOattributes.tsv')
    gvf_info_attributes_file = os.path.join(etc_folder, 'gvfINFOattributes.tsv')
    symbolic_allele_file = os.path.join(etc_folder, 'svALTkeys.tsv')

    dgva_attribute_dict = read_info_attributes(info_attributes_file=dgva_info_attributes_file) # needed to generate custom strings
    gvf_attribute_dict = read_info_attributes(info_attributes_file=gvf_info_attributes_file)
    symbolic_allele_dictionary = read_sequence_ontology_symbolic_allele(symbolic_allele_file)

    if args.assembly:
        assembly_file = os.path.abspath(args.assembly)
    else:
        assembly_file = None
    vcf_data_lines, list_of_vcf_objects = gvf_features_to_vcf_objects(gvf_lines_obj_list,
                                                                      dgva_attribute_dict,
                                                                      gvf_attribute_dict,
                                                                      symbolic_allele_dictionary,
                                                                      assembly_file,
                                                                      lines_custom_structured,
                                                                      standard_lines_dictionary,
                                                                      all_possible_lines_dictionary)


    print("Writing to the following VCF output: ", args.vcf_output)
    print("Generating the VCF header and the meta-information lines")

    with open(args.vcf_output, "w") as vcf_output:
        unique_pragmas_to_add, samples, unique_alt_lines_to_add, unique_info_lines_to_add, unique_filter_lines_to_add, unique_format_lines_to_add = generate_vcf_metainformation(lines_custom_unstructured, gvf_pragmas, gvf_non_essential, list_of_vcf_objects, standard_lines_dictionary)
        for pragma in unique_pragmas_to_add:
            vcf_output.write(f"{pragma}\n")
        for alt_lines in unique_alt_lines_to_add:
            vcf_output.write(f"{alt_lines}\n")
        for info_lines in unique_info_lines_to_add:
            vcf_output.write(f"{info_lines}\n")
        for filter_lines in unique_filter_lines_to_add:
            vcf_output.write(f"{filter_lines}\n")
        for format_lines in unique_format_lines_to_add:
            vcf_output.write(f"{format_lines}\n")
        header_fields = generate_vcf_header_line(samples)
        vcf_output.write(f"{header_fields}\n")
        print("Generating the VCF datalines")
        formatted_vcf_datalines = format_vcf_datalines(list_of_vcf_objects, samples)
        for line in formatted_vcf_datalines:
            vcf_output.write(f"{line}\n")
    vcf_output.close()
    print("GVF to VCF conversion complete")


if __name__ == "__main__":
    main()
