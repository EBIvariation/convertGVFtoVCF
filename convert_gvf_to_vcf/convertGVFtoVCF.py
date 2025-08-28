import argparse


# step 3
def generate_custom_structured_metainfomation_line(lines_custom_structured,
                                                   vcfkey,
                                                   vcfkey_id,
                                                   vcfkey_number,
                                                   vcfkey_type,
                                                   vcfkey_description,
                                                   optional_extrafields=None):
    """ Generates a custom structured meta-information line for INFO/FILTER/FORMAT/ALT
    :param vcfkey: required field INFO, FILTER, FORMAT, ALT
    :param vcfkey_id: required field for structured lines ID
    :param vcfkey_number: The number of values that can be included or special character: A or R or G or .
    :param vcfkey_type: Values are Integer, Float, Character, String
    :param vcfkey_description: Description
    :param lines_custom_structured:
    :param optional_extrafields: an optional field, dictionary of custom fields and their values
    :return: custom_structured_string
    """
    extrakeys_kvlines = []
    if optional_extrafields:
        for extra_field in optional_extrafields:
            kv_line = "," + extra_field + "=" + "\"" + optional_extrafields[extra_field] + "\""
            extrakeys_kvlines.append(kv_line)
    vcfkey_extrakeys = ''.join(extrakeys_kvlines)
    custom_structured_string = f"##{vcfkey}=<ID=\"{vcfkey_id}\",Number=\"{vcfkey_number}\",Type=\"{vcfkey_type}\",Description=\"{vcfkey_description}\"{vcfkey_extrakeys}>"
    lines_custom_structured.append(custom_structured_string)
    return custom_structured_string

# extra functions for step 4
# for INFO
def read_reserved_info_key(all_possible_INFO_lines, reservedinfokeysfile="convert_gvf_to_vcf/etc/ReservedINFOkeys.txt"):
    """ Reads in the reserved INFO keys and returns a list of all_possible_INFO_lines which can be used to populate the header

    :param reservedinfokeysfile: File to a tab-delimited table of Reserved INFO keys in Table 1 of VCF specification
    :return: all_possible_INFO_lines
    """

    with open(reservedinfokeysfile) as infokeysfile:
        next(infokeysfile)
        info_keys_content = infokeysfile.readlines()
        for info in info_keys_content:
            info_tokens = info.rstrip().split("\t")
            keyid = info_tokens[0]
            number=info_tokens[1]
            type_for_key=info_tokens[2]
            description=info_tokens[3]
            reserved_info_string = f"##INFO=<ID=\"{keyid}\",Number=\"{number}\",Type=\"{type_for_key}\",Description=\"{description}\">"
            #all_possible_INFO_lines.append(reserved_info_string)
            all_possible_INFO_lines[keyid] = reserved_info_string
    return all_possible_INFO_lines

def read_sv_info_key(all_possible_INFO_lines, svinfokeysfile="svINFOkeys.txt"):
    """ Reads in INFO keys for structural variants and return a list of all_possible_INFO_lines

    :param svinfokeysfile: File to tab delimited table of INFO keys used in Structural Variants and their VCF header
    :return: all_possible_INFO_lines
    """
    with open(svinfokeysfile) as svinfokeys:
        next(svinfokeys)
        sv_info_keys_content = svinfokeys.readlines()
        for svinfo in sv_info_keys_content:
            svinfo_tokens = svinfo.rstrip().split("\t")
            svkeyid = svinfo_tokens[0]
            svinfoline = svinfo_tokens[1]
            all_possible_INFO_lines[svkeyid]= svinfoline
    return all_possible_INFO_lines

# for FORMAT
def read_reserved_format_key(all_possible_FORMAT_lines, reservedformatkeysfile="ReservedFORMATkeys.txt"):
    """ Reads in the reserved FORMAT keys and returns a list of all_possible_FORMAT_lines which can be used to populate the header

    :param reservedformatkeysfile: file that is a tab delimited table of reserved FORMAT keys in Table 2 of VCF specification
    :return:
    """
    with open(reservedformatkeysfile) as formatkeysfile:
        next(formatkeysfile)
        format_keys_content = formatkeysfile.readlines()
        for format_key_line in format_keys_content:
            format_tokens = format_key_line.rstrip().split("\t")
            keyid = format_tokens[0]
            number = format_tokens[1]
            type_for_key = format_tokens[2]
            description = format_tokens[3]
            reserved_format_string = f"##FORMAT=<ID=\"{keyid}\",Number=\"{number}\",Type=\"{type_for_key}\",Description=\"{description}\">"
            all_possible_FORMAT_lines[keyid] = reserved_format_string
    return all_possible_FORMAT_lines

def read_sv_format_keys(all_possible_FORMAT_lines, svformatkeysfile="svFORMATkeys.txt"):
    """ Reads in FORMAT keys for strucural variants and returns a list of all_possible_FORMAT_lines

    :param svformatkeysfile: File to tab delimited table of FORMAT keys used in Structural Variants and their VCF header
    :return: all_possible_FORMAT_lines
    """
    with open(svformatkeysfile) as svformatkeys:
        next(svformatkeys)
        sv_format_keys_content = svformatkeys.readlines()
        for svformat in sv_format_keys_content:
            svformat_tokens = svformat.rstrip().split()
            svkeyid = svformat_tokens[0]
            svformatline = svformat_tokens[1]
            all_possible_FORMAT_lines[svkeyid] = svformatline
    return all_possible_FORMAT_lines

# for ALT
def read_sv_alt_keys(all_possible_ALT_lines, svaltkeysfile="svALTkeys.txt"):
    """ Reads in ALT keys for structural variants and return a list of all_possible_ALT_lines

    :param svaltkeysfile: File to tab delimited table of ALT keys used in Structural Variants and their VCF header
    :return: all_possible_ALT_lines
    """
    with open(svaltkeysfile) as svaltkeys:
        next(svaltkeys)
        sv_alt_keys_content = svaltkeys.readlines()
        for svalt in sv_alt_keys_content:
            svalt_tokens = svalt.rstrip().split()
            svkeyid = svalt_tokens[0]
            svaltline = svalt_tokens[1]
            all_possible_ALT_lines[svkeyid] = svaltline
        return all_possible_ALT_lines

def generate_all_standard_structured_metainformation_line(vcfkey, all_possible_ALT_lines, all_possible_INFO_lines, all_possible_FILTER_lines, all_possible_FORMAT_lines):
    #, vcfkey_id,vcf_number, vcf_type, vcf_description):
    if vcfkey=="INFO":
        # generate all possible lines for the reserved info keys
        read_reserved_info_key(all_possible_INFO_lines, reservedinfokeysfile="convert_gvf_to_vcf/etc/ReservedINFOkeys.txt")
        read_sv_info_key(all_possible_INFO_lines, svinfokeysfile="convert_gvf_to_vcf/etc/svINFOkeys.txt")
        return all_possible_ALT_lines, all_possible_INFO_lines, all_possible_FILTER_lines, all_possible_FORMAT_lines
    elif vcfkey=="FORMAT":
        # TABLE 2
        # FORMAT KEYS FOR STRUCTURAL VARIANTS
        read_reserved_format_key(all_possible_FORMAT_lines, reservedformatkeysfile="convert_gvf_to_vcf/etc/ReservedFORMATkeys.txt")
        read_sv_format_keys(all_possible_FORMAT_lines, svformatkeysfile="convert_gvf_to_vcf/etc/svFORMATkeys.txt")
        return all_possible_ALT_lines, all_possible_INFO_lines, all_possible_FILTER_lines, all_possible_FORMAT_lines
    elif vcfkey=="FILTER":
        # MAY NOT BE NEEDED
        pass
    elif vcfkey=="ALT":
        # note: svALTkey may be an incomplete list at the moment
        # no reserved alt keys
        read_sv_alt_keys(all_possible_ALT_lines, svaltkeysfile="convert_gvf_to_vcf/etc/svALTkeys.txt")
        return all_possible_ALT_lines, all_possible_INFO_lines, all_possible_FILTER_lines, all_possible_FORMAT_lines
    else:
        print("Please provide a key: INFO, FORMAT,FILTER, ALT")
        return None

# step 4
def generate_standard_structured_metainformation_line(vcfkey, vcfkeyid, lines_standard_ALT, lines_standard_INFO, lines_standard_FILTER, lines_standard_FORMAT, all_possible_ALT_lines, all_possible_INFO_lines, all_possible_FILTER_lines, all_possible_FORMAT_lines):
    standard_structured_line = ""
    if vcfkey=="INFO":
        # print("# A start")
        # print(all_possible_INFO_lines)
        # print("# A end")
        if vcfkeyid in all_possible_INFO_lines:
            standard_structured_line = all_possible_INFO_lines[vcfkeyid]
            lines_standard_INFO.append(standard_structured_line)
            # print(standard_structured_line)
    elif vcfkey=="FORMAT":
        standard_structured_line = all_possible_FORMAT_lines[vcfkeyid]
        lines_standard_FORMAT.append(standard_structured_line)
    elif vcfkey=="FILTER":
        standard_structured_line = all_possible_FILTER_lines[vcfkeyid]
        lines_standard_FILTER.append(standard_structured_line)
    elif vcfkey=="ALT":
        standard_structured_line = all_possible_ALT_lines[vcfkeyid]
        lines_standard_ALT.append(standard_structured_line)
    else:
        print("Please provide a key: INFO, FORMAT,FILTER, ALT")
    return standard_structured_line

# step 5
def generate_custom_unstructured_metainfomation_line(vcf_unstructured_key, vcf_unstructured_value, lines_custom_unstructured):
    """ Generates a formatted unstructured metainformation line using a custom key value pair. This is stored in the list called lines_custom_unstructured.
    :param lines_custom_unstructured: list to store custom unstructured metainformation lines
    :param vcf_unstructured_key: key for custom unstructured metainformation line
    :param vcf_unstructured_value: value for custom unstructured metainformation line
    :return: custom_unstructured_string
    """
    custom_unstructured_string = f"##{vcf_unstructured_key}={vcf_unstructured_value}"
    lines_custom_unstructured.append(custom_unstructured_string)
    return custom_unstructured_string

# additional support for writing to file
def write_to_vcf_file(vcf_output_file, string_to_add):
    with open(vcf_output_file, "a") as vcf:
        vcf.write(string_to_add)
        vcf.write("\n")

# additional support function for step 6
def read_dgva_info_attributes(dgva_info_attributes_file):
    """ Read in the file containing DGVa specific INFO attributes.
    :param dgva_info_attributes_file: A file containing the DGVa specific attributes
    :return: dgva_attribute_dict: A dictionary of key id and the list of attribute tokens
    """
    dgva_attribute_dict = {}  # dictionary of dgva specific INFO attributes
    with open(dgva_info_attributes_file) as dgva_info_attributes:
        next(dgva_info_attributes)
        dgva_info_attributes_contents = dgva_info_attributes.readlines()
        for dgva_attribute in dgva_info_attributes_contents:
            dgva_attribute_tokens = dgva_attribute.rstrip().split("\t")
            dgva_key_id = dgva_attribute_tokens[0]
            dgva_attribute_dict[dgva_key_id] = dgva_attribute_tokens
    return dgva_attribute_dict

# additional support function for step 6
def read_gvf_info_attributes(gvf_info_attributes_file):
    """ Read in the file of GVF INFO attributes
    :param gvf_info_attributes_file: file of GVF INFO attributes
    :return: gvf_attribute_dict: a dictionary of gvf attributes
    """
    gvf_attribute_dict = {}  # dictionary of GVF attributes
    with open(gvf_info_attributes_file) as gvf_info_attributes:
        next(gvf_info_attributes)
        gvf_info_attributes_contents = gvf_info_attributes.readlines()
        for gvf_attribute in gvf_info_attributes_contents:
            gvf_attribute_tokens = gvf_attribute.rstrip().split("\t")
            gvf_key_id = gvf_attribute_tokens[0]
            gvf_attribute_dict[gvf_key_id] = gvf_attribute_tokens
    return gvf_attribute_dict

def get_gvf_attributes(column9_of_gvf):
    """Get a dictionary of GVF attributes
    :param column9_of_gvf: attributes column - the final column of the GVF file
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
                                         lines_standard_ALT,
                                         lines_standard_INFO,
                                         lines_standard_FILTER,
                                         lines_standard_FORMAT,
                                         all_possible_ALT_lines,
                                         all_possible_INFO_lines,
                                         all_possible_FILTER_lines,
                                         all_possible_FORMAT_lines):

    gvf_attribute_dictionary = get_gvf_attributes(column9_of_gvf)
    vcf_vals = {}
    catching_for_review = []
    # created a rough guide to attributes_for_custom_structured_metainfomation in dgvaINFOattributes.txt = this probably should be refined at a later date
    # TODO: edit dgvaINFOattributes.txt i.e. replace unknown placeholders '.' with the actual answer, provide a more informative description
    for attrib_key in gvf_attribute_dictionary:
        # if dgva specific key, create custom string otherwise do standard
        if attrib_key in dgva_attribute_dict:
            generate_custom_structured_metainfomation_line(lines_custom_structured,
                                                           vcfkey="INFO", vcfkey_id=attrib_key,
                                                           vcfkey_number=dgva_attribute_dict[attrib_key][1],
                                                           vcfkey_type=dgva_attribute_dict[attrib_key][2],
                                                           vcfkey_description=dgva_attribute_dict[attrib_key][3],
                                                           optional_extrafields=None)
            vcf_vals[attrib_key]=gvf_attribute_dictionary[attrib_key]
        elif attrib_key == "allele_count":
            generate_standard_structured_metainformation_line("INFO", "AC", lines_standard_ALT, lines_standard_INFO, lines_standard_FILTER, lines_standard_FORMAT, all_possible_ALT_lines, all_possible_INFO_lines, all_possible_FILTER_lines, all_possible_FORMAT_lines)
            vcf_vals[attrib_key]=gvf_attribute_dictionary[attrib_key]
            # print("# VCF start")
            # print([-1])
            # print(vcf_vals["allele_count"])
            # print("# VCF END")
        elif attrib_key == "allele_frequency":
            generate_standard_structured_metainformation_line("INFO", "AF", lines_standard_ALT, lines_standard_INFO, lines_standard_FILTER, lines_standard_FORMAT, all_possible_ALT_lines, all_possible_INFO_lines, all_possible_FILTER_lines, all_possible_FORMAT_lines)
            # print(attrib_key + "=" + str(gvf_attribute_dictionary[attrib_key])) # add this to the INFO
        elif attrib_key == "ciend":
            generate_standard_structured_metainformation_line("INFO", "CIEND", lines_standard_ALT, lines_standard_INFO, lines_standard_FILTER, lines_standard_FORMAT, all_possible_ALT_lines, all_possible_INFO_lines, all_possible_FILTER_lines, all_possible_FORMAT_lines)
            # print(attrib_key + "=" + str(gvf_attribute_dictionary[attrib_key])) # add this to the INFO
        elif attrib_key == "copy_number":
            generate_standard_structured_metainformation_line("INFO", "CN", lines_standard_ALT, lines_standard_INFO, lines_standard_FILTER, lines_standard_FORMAT, all_possible_ALT_lines, all_possible_INFO_lines, all_possible_FILTER_lines, all_possible_FORMAT_lines)
            # print(attrib_key + "=" + str(gvf_attribute_dictionary[attrib_key])) # add this to the INFO
        elif attrib_key == "insertion_length":
            generate_standard_structured_metainformation_line("INFO", "SVLEN", lines_standard_ALT, lines_standard_INFO, lines_standard_FILTER, lines_standard_FORMAT, all_possible_ALT_lines, all_possible_INFO_lines, all_possible_FILTER_lines, all_possible_FORMAT_lines)
            # print(attrib_key + "=" + str(gvf_attribute_dictionary[attrib_key])) # add this to the INFO
        elif attrib_key == "mate_id":
            generate_standard_structured_metainformation_line("INFO","MATEID", lines_standard_ALT, lines_standard_INFO, lines_standard_FILTER, lines_standard_FORMAT, all_possible_ALT_lines, all_possible_INFO_lines, all_possible_FILTER_lines, all_possible_FORMAT_lines)
            # print(attrib_key + "=" + str(gvf_attribute_dictionary[attrib_key])) # add this to the INFO
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
            generate_custom_structured_metainfomation_line(lines_custom_structured,vcfkey="INFO", vcfkey_id=attrib_key,
                                                           vcfkey_number=gvf_attribute_dict[attrib_key][1],
                                                           vcfkey_type=gvf_attribute_dict[attrib_key][2],
                                                           vcfkey_description=gvf_attribute_dict[attrib_key][3],
                                                           optional_extrafields=None)
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
        :param source: source i.e DGVa
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

    def print_gvffeatureline(self):
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
class VcfDataObj:
    def __init__(self, gvf_feature_line_object,
                 dgva_attribute_dict,
                 gvf_attribute_dict,
                 lines_custom_structured,
                 lines_standard_ALT,
                 lines_standard_INFO,
                 lines_standard_FILTER,
                 lines_standard_FORMAT,
                 all_possible_ALT_lines,
                 all_possible_INFO_lines,
                 all_possible_FILTER_lines,
                 all_possible_FORMAT_lines):
        # ATTRIBUTES
        self.vcf_value = convert_gvf_attributes_to_vcf_values(gvf_feature_line_object.attributes,
                                                              dgva_attribute_dict,
                                                              gvf_attribute_dict,
                                                              lines_custom_structured,
                                                              lines_standard_ALT,
                                                              lines_standard_INFO,
                                                              lines_standard_FILTER,
                                                              lines_standard_FORMAT,
                                                              all_possible_ALT_lines,
                                                              all_possible_INFO_lines,
                                                              all_possible_FILTER_lines,
                                                              all_possible_FORMAT_lines)

        # DATALINE
        self.chrom = gvf_feature_line_object.seqid
        self.pos = gvf_feature_line_object.start
        self.id = self.vcf_value["ID"]  # attributes: ID
        self.ref = self.get_ref()
        self.alt = self.vcf_value["Variant_seq"] # attributes: variant_seq
        self.qual = gvf_feature_line_object.score # see EVA-3879: this is always '.'
        self.filter = "." # this is always a placeholder '.'; perhaps could add s50.

        # INFO
        #TODO: specific SV info keys populated from gvf_feature_line
        self.key = self.chrom + "_" + self.pos
        self.info = "pending_aggregation" # TODO: add info field for self.info
        self.source = gvf_feature_line_object.source

        # # non-vcf
        self.phase = gvf_feature_line_object.phase # this is always a placeholder'.'
        self.end = gvf_feature_line_object.end
        self.so_type = gvf_feature_line_object.feature_type #currently column 3 of gvf, but could be an attribute so perhapsVCF: INFO or FORMAT?

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

    def get_ref(self):
        """ Gets the reference allele from attributes column or if not found, returns "."
        :return: reference allele
        """
        if "Reference_seq" in self.vcf_value.keys():
            return self.vcf_value["Reference_seq"] # attributes:reference_seq
        else:
            return "." # TODO: how shall we fill this in with this scenario?

    def get_string(self):
        string_to_return = (self.chrom + "\t" +
                            self.pos + "\t" +
                            self.key + "\t" +
                            self.qual + "\t" +
                            self.filter + "\t" +
                            self.info + "\t" +
                            self.source + "\t" +
                            self.phase + "\t" +
                            self.end + "\t" +
                            self.so_type + "\t" +
                            self.sample_name + "\t" +
                            self.format
                            )
        return string_to_return

#step 9 using custom unstructured meta-information line = generate_custom_unstructured_metainfomation_line
def generate_vcf_metainformation(lines_custom_unstructured, gvf_pragmas, list_of_vcf_objects):
    """ Generates a list of metainformation lines for the VCF header
    :param lines_custom_unstructured: a list of formatted unstructured metainformation lines using a custom key value pair
    :param gvf_pragmas: list of gvf pragmas to convert
    :param list_of_vcf_objects: list of vcf objects
    :return: unique_pragmas_to_add: a list of pragmas, this list contains no duplicates
    """
    pragmas_to_add = []
    unique_pragmas_to_add = []
    # MANDATORY: file format for VCF
    pragma_fileformat = generate_custom_unstructured_metainfomation_line("fileformat", "VCFv4.4",lines_custom_unstructured)
    pragmas_to_add.append(pragma_fileformat)
    # print("# PRAGMA start 1")
    # print(pragma_fileformat)
    # print(pragmas_to_add[-1])
    # print("# PRAGMA end 1")
    #TODO: list of pragmas to add:reference=file, contig, phasing,INFO#
    for pragma in gvf_pragmas:
        # file date
        if pragma.startswith("##file-date"):
            date = pragma.split(" ")[1].replace("-", "")
            pragma_filedate = generate_custom_unstructured_metainfomation_line("fileDate", date, lines_custom_unstructured)
            pragmas_to_add.append(pragma_filedate)
        # source
        for vcf_obj in list_of_vcf_objects:
            pragma_source = generate_custom_unstructured_metainfomation_line("source", vcf_obj.source, lines_custom_unstructured)
            pragmas_to_add.append(pragma_source)
        # reference #TODO: add this
        # contig: recommended, see section 1.4.7 of VCF specification # TODO: add this
        # phasing # not required
        if pragma.startswith("##gff-version"):
            gff_version_number = pragma.split(" ")[1]
            pragma_gff_version_number = generate_custom_unstructured_metainfomation_line("gff-version", gff_version_number, lines_custom_unstructured)
            pragmas_to_add.append(pragma_gff_version_number)
        elif pragma.startswith("##gvf-version"):
            gvf_version_number = pragma.split(" ")[1]
            pragma_gvf_version_number = generate_custom_unstructured_metainfomation_line("gvf-version", gvf_version_number, lines_custom_unstructured)
            pragmas_to_add.append(pragma_gvf_version_number)
        elif pragma.startswith("##species"):
            species_value = pragma.split(" ")[1]
            pragma_species_value = generate_custom_unstructured_metainfomation_line("species", species_value, lines_custom_unstructured)
            pragmas_to_add.append(pragma_species_value)
        elif pragma.startswith("##genome-build"):
            genome_build = pragma.split("genome-build ")[1]
            pragma_genome_build = generate_custom_unstructured_metainfomation_line("genome-build", genome_build, lines_custom_unstructured)
            pragmas_to_add.append(pragma_genome_build)
        else:
            pass
    for pragma in pragmas_to_add:
        if pragma not in unique_pragmas_to_add:
            unique_pragmas_to_add.append(pragma)

    #TODO: add the pragmas for the GVF non-essentials
    return unique_pragmas_to_add

# step 10
# TODO: finish the below for sample names
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
                                lines_custom_structured,
                                lines_standard_ALT,
                                lines_standard_INFO,
                                lines_standard_FILTER,
                                lines_standard_FORMAT,
                                all_possible_ALT_lines,
                                all_possible_INFO_lines,
                                all_possible_FILTER_lines,
                                all_possible_FORMAT_lines):
    """ Creates VCF objects from GVF feature lines and stores the VCF objects.
    :param gvf_lines_obj_list: list of GVF feature line objects
    :param gvf_attribute_dict: dictionary of GVF INFO attributes
    :param dgva_attribute_dict: dictionary af DGVa specific INFO attributes
    :param lines_custom_structured: list to store custom structured metainformation lines
    :param lines_standard_ALT: ALT lines for this VCF file
    :param lines_standard_INFO: INFO lines for this VCF file
    :param lines_standard_FILTER: FILTER lines for this VCF file
    :param lines_standard_FORMAT: FORMAT lines for this VCF file
    :param all_possible_ALT_lines: dict of all possible ALT lines
    :param all_possible_INFO_lines: dict of all possible INFO lines
    :param all_possible_FILTER_lines: dict of all possible FILTER lines
    :param all_possible_FORMAT_lines: dict of all possible FORMAT lines
    :return: vcf_data_lines, list_of_vcf_objects: dictionary of lists and a list of VCF objects
    """
    vcf_data_lines = {}  # DICTIONARY OF LISTS
    list_of_vcf_objects = []

    # create a vcf object for every feature line in the GVF (1:1)
    # add the newly created vcf object to the vcf data line it belongs to
    # (1:many; key=chrom_pos; 1 key: many vcf objects)
    for gvf_featureline in gvf_lines_obj_list:
        vcf_object = VcfDataObj(gvf_featureline,
                                dgva_attribute_dict,
                                gvf_attribute_dict,
                                lines_custom_structured,
                                lines_standard_ALT,
                                lines_standard_INFO,
                                lines_standard_FILTER,
                                lines_standard_FORMAT,
                                all_possible_ALT_lines,
                                all_possible_INFO_lines,
                                all_possible_FILTER_lines,
                                all_possible_FORMAT_lines)
        list_of_vcf_objects.append(vcf_object)
        if vcf_object.key in vcf_data_lines:
            vcf_data_lines[vcf_object.key].append(vcf_object)
        else:
            vcf_data_line_objects_list = []
            vcf_data_line_objects_list.append(vcf_object)
            vcf_data_lines[vcf_object.key] = vcf_data_line_objects_list
        # print("# start - creating vcf object - 1")
        # print(vcf_data_lines[vcf_object.key][-1].get_string())
        # print("# end - creating vcf object - 1")
        # print("number of keys in dict", len(vcf_data_lines.keys()))
        # check the number of objects to see if they are merged
        # for key in vcf_data_lines.keys():
        #     vcf_obj_list = vcf_data_lines[key]
            # print("for ", key, " the number of vcf objects is: ", len(vcf_obj_list))
    return vcf_data_lines, list_of_vcf_objects

def format_vcf_datalines(list_of_vcf_objects):
    """ Iterates through a list of VCF objects and formats them as a VCF dataline.
    :param list_of_vcf_objects: list of vcf objects
    :return: formatted_vcf_datalines: list of formatted vcf datalines
    """
    formatted_vcf_datalines = []
    for vcf_obj in list_of_vcf_objects:
        vcf_line = (f"{vcf_obj.chrom}\t"
                    f"{vcf_obj.pos}\t"
                    f"{vcf_obj.id}\t"
                    f"{vcf_obj.ref}\t" #TODO: should this always be empty
                    f"{vcf_obj.alt}\t" #TODO: should this always be empty
                    f"{vcf_obj.qual}\t" #TODO: should this always be empty
                    f"{vcf_obj.filter}\t" #TODO: should this always be empty
                    f"{vcf_obj.info}\t"
                    f"{vcf_obj.format}\tsampleFORMAThere" #TODO: fill this in
                    )
        formatted_vcf_datalines.append(vcf_line)
    return formatted_vcf_datalines

def main():
    print("Running the GVF to VCF converter")
    # step 1
    parser = argparse.ArgumentParser()
    parser.add_argument("gvf_input", help="GVF input file.")
    parser.add_argument("vcf_output", help="VCF output file.")
    args = parser.parse_args()
    print("The provided input file is: ", args.gvf_input)
    print("The provided output file is: ", args.vcf_output)
    # Dictionary for all possible VCF meta-information lines
    all_possible_ALT_lines = {}
    all_possible_INFO_lines = {}  # dictionary, ID => INFO meta-information line for that particular ID
    all_possible_FILTER_lines = {}
    all_possible_FORMAT_lines = {}  # dictionary, ID => FORMAT meta-information line for that particular ID

    all_possible_ALT_lines, all_possible_INFO_lines, all_possible_FILTER_lines, all_possible_FORMAT_lines = generate_all_standard_structured_metainformation_line("INFO", all_possible_ALT_lines, all_possible_INFO_lines,all_possible_FILTER_lines, all_possible_FORMAT_lines)
    all_possible_ALT_lines, all_possible_INFO_lines, all_possible_FILTER_lines, all_possible_FORMAT_lines = generate_all_standard_structured_metainformation_line("ALT", all_possible_ALT_lines, all_possible_INFO_lines,all_possible_FILTER_lines, all_possible_FORMAT_lines)
    all_possible_ALT_lines, all_possible_INFO_lines, all_possible_FILTER_lines, all_possible_FORMAT_lines = generate_all_standard_structured_metainformation_line("FORMAT", all_possible_ALT_lines, all_possible_INFO_lines, all_possible_FILTER_lines, all_possible_FORMAT_lines)
    # custom meta-information lines for this VCF file
    lines_custom_structured = []
    lines_custom_unstructured = []
    # standard structured meta-information lines for this VCF file
    lines_standard_ALT = []
    lines_standard_INFO = []
    lines_standard_FILTER = []
    lines_standard_FORMAT = []


    gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(args.gvf_input)
    dgva_attribute_dict = read_dgva_info_attributes(dgva_info_attributes_file="convert_gvf_to_vcf/etc/dgvaINFOattributes.txt") # needed to generate custom strings
    gvf_attribute_dict = read_gvf_info_attributes(gvf_info_attributes_file="convert_gvf_to_vcf/etc/gvfINFOattributes.txt")
    vcf_data_lines, list_of_vcf_objects = gvf_features_to_vcf_objects(gvf_lines_obj_list,
                                                                      dgva_attribute_dict,
                                                                      gvf_attribute_dict,
                                                                      lines_custom_structured,
                                                                      lines_standard_ALT,
                                                                      lines_standard_INFO,
                                                                      lines_standard_FILTER,
                                                                      lines_standard_FORMAT,
                                                                      all_possible_ALT_lines,
                                                                      all_possible_INFO_lines,
                                                                      all_possible_FILTER_lines,
                                                                      all_possible_FORMAT_lines)

    # 10c
    print("Writing to the following VCF output: ", args.vcf_output)
    print("Generating the VCF header and the meta-information lines")
    unique_pragmas_to_add = generate_vcf_metainformation(lines_custom_unstructured, gvf_pragmas, list_of_vcf_objects)
    for pragma in unique_pragmas_to_add:
        write_to_vcf_file(args.vcf_output, pragma)
    samples = ["samA"] # TODO: this is a placeholder, need to add a function to read gvf pragmas and collect the samples into a list
    header_fields = generate_vcf_header_line(samples)
    write_to_vcf_file(args.vcf_output, header_fields)
    print("Generating the VCF datalines")
    formatted_vcf_datalines = format_vcf_datalines(list_of_vcf_objects)
    for line in formatted_vcf_datalines:
        write_to_vcf_file(args.vcf_output, line)
    print("GVF to VCF conversion complete")


if __name__ == "__main__":
    main()
