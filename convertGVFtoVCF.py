import argparse

# unforeseen extra structures
# Dictionary for all possible VCF meta-information lines
all_possible_ALT_lines = {}
all_possible_INFO_lines = {} # dictionary, ID => INFO meta-information line for that particular ID
all_possible_FILTER_lines = {}
all_possible_FORMAT_lines = {} # dictionary, ID => FORMAT meta-information line for that particular ID
# standard structured meta-information lines for this VCF file
lines_ALT = []
lines_INFO = []
lines_FILTER = []
lines_FORMAT = []
# custom meta-information lines for this VCF file
LINES_CUSTOM_STRUCTURED = [] # potential for source
lines_custom_unstructured = [] # meta-information lines like ##source

sample_names = []

# global lists for reading in the gvf file
GVF_PRAGMAS = [] # list of pragma lines starting with: ##
GVF_NON_ESSENTIAL = [] # list of non-essential lines starting with: #
GVF_LINES_OBJ_LIST = [] # list of objects when reading in gvf files, one object represents a gvf line

# global dictionary for storing attributes which are specific to DGVa
DGVA_ATTRIBUTE_DICT = {} # dictionary of dgva INFO attributes

# global dictionary for storing GVF attributes
GVF_ATTRIBUTE_DICT = {} # dictionary of GVF attributes

# step 3
def generate_custom_structured_metainfomation_line(vcfkey,
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
    LINES_CUSTOM_STRUCTURED.append(custom_structured_string)
    return custom_structured_string

# unforeseen extra functions for step 4
# for INFO
def read_reserved_info_key(reservedinfokeysfile="ReservedINFOkeys.txt"):
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

def read_sv_info_key(svinfokeysfile="svINFOkeys.txt"):
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
def read_reserved_format_key(reservedformatkeysfile="ReservedFORMATkeys.txt"):
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
def read_sv_format_keys(svformatkeysfile="svFORMATkeys.txt"):
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
def read_sv_alt_keys(svaltkeysfile="svALTkeys.txt"):
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

def generate_all_standard_structured_metainformation_line(vcfkey):
    #, vcfkey_id,vcf_number, vcf_type, vcf_description):
    if vcfkey=="INFO":
        # generate all possible lines for the reserved info keys
        read_reserved_info_key(reservedinfokeysfile="ReservedINFOkeys.txt")
        read_sv_info_key(svinfokeysfile="svINFOkeys.txt")
        return all_possible_INFO_lines
    elif vcfkey=="FORMAT":
        # TABLE 2
        # FORMAT KEYS FOR STRUCTURAL VARIANTS
        read_reserved_format_key(reservedformatkeysfile="ReservedFORMATkeys.txt")
        read_sv_format_keys(svformatkeysfile="svFORMATkeys.txt")
        return all_possible_FORMAT_lines
    elif vcfkey=="FILTER":
        # MAY NOT BE NEEDED
        pass
    elif vcfkey=="ALT":
        # note: svALTkey may be an incomplete list at the moment
        # no reserved alt keys
        read_sv_alt_keys(svaltkeysfile="svALTkeys.txt")
        return all_possible_ALT_lines
    else:
        print("Please provide a key: INFO, FORMAT,FILTER, ALT")
        return None
# step 4
def generate_standard_structured_metainformation_line(vcfkey, vcfkeyid):
    standard_structured_line = ""
    if vcfkey=="INFO":
        standard_structured_line = all_possible_INFO_lines[vcfkeyid]
        lines_INFO.append(standard_structured_line)
        print(standard_structured_line)
    elif vcfkey=="FORMAT":
        standard_structured_line = all_possible_FORMAT_lines[vcfkeyid]
        lines_FORMAT.append(standard_structured_line)
    elif vcfkey=="FILTER":
        standard_structured_line = all_possible_FILTER_lines[vcfkeyid]
        lines_FILTER.append(standard_structured_line)
    elif vcfkey=="ALT":
        standard_structured_line = all_possible_ALT_lines[vcfkeyid]
        lines_ALT.append(standard_structured_line)
    else:
        print("Please provide a key: INFO, FORMAT,FILTER, ALT")
    return standard_structured_line

# step 5
def generate_custom_unstructured_metainfomation_line(vcf_unstructured_key, vcf_unstructured_value):
    """ Generates a formatted unstructured metainformation line using a custom key value pair. This is stored in the list called lines_custom_unstructured.
    :param vcf_unstructured_key: key for custom unstructured metainformation line
    :param vcf_unstructured_value: value for custom unstructured metainformation line
    :return: custom_unstructured_string
    """
    custom_unstructured_string = f"##{vcf_unstructured_key}={vcf_unstructured_value}"
    lines_custom_unstructured.append(custom_unstructured_string)
    return custom_unstructured_string

# additional support function for step 6
def read_dgva_info_attributes(dgva_info_attributes_file="dgvaINFOattributes.txt"):
    """ Read in the file containing DGVa specific INFO attributes.
    :param dgva_info_attributes_file: A file containing the DGVa specific attributes
    :return: None
    """
    with open(dgva_info_attributes_file) as dgva_info_attributes:
        next(dgva_info_attributes)
        dgva_info_attributes_contents = dgva_info_attributes.readlines()
        for dgva_attribute in dgva_info_attributes_contents:
            dgva_attribute_tokens = dgva_attribute.rstrip().split("\t")
            dgva_key_id = dgva_attribute_tokens[0]
            DGVA_ATTRIBUTE_DICT[dgva_key_id] = dgva_attribute_tokens


# additional support function for step 6
def read_gvf_info_attributes(gvf_info_attributes_file="gvfINFOattributes.txt"):
    """ Read in the file of GVF INFO attributes
    :param gvf_info_attributes_file: file of GVF INFO attributes
    :return: None
    """
    with open(gvf_info_attributes_file) as gvf_info_attributes:
        next(gvf_info_attributes)
        gvf_info_attributes_contents = gvf_info_attributes.readlines()
        for gvf_attribute in gvf_info_attributes_contents:
            gvf_attribute_tokens = gvf_attribute.rstrip().split("\t")
            gvf_key_id = gvf_attribute_tokens[0]
            GVF_ATTRIBUTE_DICT[gvf_key_id] = gvf_attribute_tokens



catching_for_review = []
# step 6
# CAVEATS: 1) assume sample_name is present in the GVF file. If absent consider adding UnknownSample1, UnknownSample2 etc.
def convert_gvf_attributes(column9_of_gvf):
    # PART 1: PARSE
    dictionary_of_attribute_in_gvf_line = {} # attibute key => value
    # parse by semicolon this creates attribute
    # parse by equals sign this creates tag-values, if the value is a comma, create a list
    attributes_in_gvf_line = column9_of_gvf.split(";")
    for attribute in attributes_in_gvf_line:
        attribute_key, attribute_value = attribute.split("=")
        if "," in attribute_value:
            attribute_value_list = attribute_value.split(",")
            dictionary_of_attribute_in_gvf_line[attribute_key] = attribute_value_list
        else:
            dictionary_of_attribute_in_gvf_line[attribute_key] = attribute_value
    # PART 2: loop through attributes and make a decision
    # create a custom INFO tag for all in the list attributes_for_custom_structured_metainfomation
    # below is a list of special terms for generating custom structured lines i.e. DGVa specific terminology
    attributes_for_custom_structured_metainfomation = ["3'_inner_flank_link", "5'_inner_flank_link", "allele_number", "assertion_method", "breakpoint_order",
        "clinical_significance", "evidence_sequence_link", "log2_value", "mutation_order", "Name", "parent",
        "phenotype_description", "phenotype_link", "reciprocal_alignment", "remap_score", "samples",
        "submitter_variant_call_id", "submitter_variant_region_id", "support_count", "supporting_sequence_link",
        "variant_call_description"                                                      ]
    # created a rough guide to attributes_for_custom_structured_metainfomation in dgvaINFOattributes.txt = this probably should be refined at a later date
    # TODO: edit dgvaINFOattributes.txt i.e. replace unknowns '.' with the actual answer, provide a more informative description
    for attrib_key in dictionary_of_attribute_in_gvf_line:
        # if dgva specific key, create custom string otherwise do standard
        if attrib_key in attributes_for_custom_structured_metainfomation:
            generate_custom_structured_metainfomation_line(vcfkey="INFO", vcfkey_id=attrib_key,
                                                           vcfkey_number=DGVA_ATTRIBUTE_DICT[attrib_key][1],
                                                           vcfkey_type=DGVA_ATTRIBUTE_DICT[attrib_key][2],
                                                           vcfkey_description=DGVA_ATTRIBUTE_DICT[attrib_key][3],
                                                           optional_extrafields=None)
             # add this to the INFO
        elif attrib_key == "allele_count":
            generate_standard_structured_metainformation_line("INFO", "AC")
            #print(attrib_key + "=" + str(dictionary_of_attribute_in_gvf_line[attrib_key])) # add this to the info
        elif attrib_key == "allele_frequency":
            generate_standard_structured_metainformation_line("INFO", "AF")
            # print(attrib_key + "=" + str(dictionary_of_attribute_in_gvf_line[attrib_key])) # add this to the INFO
        elif attrib_key == "ciend":
            generate_standard_structured_metainformation_line("INFO", "CIEND")
            # print(attrib_key + "=" + str(dictionary_of_attribute_in_gvf_line[attrib_key])) # add this to the INFO
        elif attrib_key == "copy_number":
            generate_standard_structured_metainformation_line("INFO", "CN")
            # print(attrib_key + "=" + str(dictionary_of_attribute_in_gvf_line[attrib_key])) # add this to the INFO
        elif attrib_key == "insertion_length":
            generate_standard_structured_metainformation_line("INFO", "SVLEN")
            # print(attrib_key + "=" + str(dictionary_of_attribute_in_gvf_line[attrib_key])) # add this to the INFO
        elif attrib_key == "mate_id":
            generate_standard_structured_metainformation_line("INFO","MATEID")
            # print(attrib_key + "=" + str(dictionary_of_attribute_in_gvf_line[attrib_key])) # add this to the INFO
        elif attrib_key == "sample_name":
            sample_names.append(sample_names)
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
            generate_custom_structured_metainfomation_line(vcfkey="INFO", vcfkey_id=attrib_key,
                                                           vcfkey_number=GVF_ATTRIBUTE_DICT[attrib_key][1],
                                                           vcfkey_type=GVF_ATTRIBUTE_DICT[attrib_key][2],
                                                           vcfkey_description=GVF_ATTRIBUTE_DICT[attrib_key][3],
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
    print("dictionary", dictionary_of_attribute_in_gvf_line)
    return dictionary_of_attribute_in_gvf_line

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
def read_in_gvf_file(arguments):
    """ Reads in the user provided GVF file.
    :param arguments: arguments.gvf_input
    :return: None
    """
    features = []
    print("Reading in the following GVF input: " + arguments.gvf_input)
    with open(arguments.gvf_input) as gvf_file:
        gvf_content = gvf_file.readlines()
    for line in gvf_content:
        if line.startswith("##"):
            GVF_PRAGMAS.append(line.rstrip())
        elif line.startswith("#"):
            GVF_NON_ESSENTIAL.append(line.rstrip())
        else:
            features.append(line.rstrip())
    for feature in features:
        f_list = feature.split("\t")
        line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6], f_list[7], f_list[8])
        GVF_LINES_OBJ_LIST.append(line_object)





#step 8
#TODO: ID this can be a semi-colon separated list or a '.' (if no value = '.'; one value = value; more than one = value;value)
class VcfDataObj:
    def __init__(self, gvf_featureline):
        self.chrom = gvf_featureline.seqid
        self.pos = gvf_featureline.start
        self.key = self.chrom + "_" + self.pos
        self.qual = gvf_featureline.score # see EVA-3879: this is always '.'
        self.filter = "." # this is always a placeholder '.'; perhaps could add s50.
        self.info = "pending aggregation"
        self.source = gvf_featureline.source
        # # non-vcf
        self.phase = gvf_featureline.phase # this is always a placeholder'.'
        self.end = gvf_featureline.end
        self.so_type = gvf_featureline.feature_type #currently column 3 of gvf, but could be an attribute so perhapsVCF: INFO or FORMAT?
        self.attributes = convert_gvf_attributes(gvf_featureline.attributes)
        self.sample_name = self.attributes["sample_name"] # this should be each samples names format value # sample names needs to be populated in attributes
        # # higher priority
        self.format = "pending" #TODO: set this in convertgvfattributes
        self.id =  self.attributes["ID"] # attributes: ID
        # if "Reference_seq" in self.attributes.keys():
        #     self.ref = self.attributes["Reference_seq"] # attributes:reference_seq
        # else:
        #     self.ref = "." # TODO: how shall we fill this in with this scenario?
        # self.alt = self.attributes["Variant_seq"]# attributes: variant_seq
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

# TODO: upto step 8, now must find a way to add info field for self.info

#step 9 using custom unstructured meta-information line = generate_custom_unstructured_metainfomation_line
def generate_vcf_metainformation():
    pragmas_to_add = []
    unique_pragmas_to_add = []
    # MANDATORY: file format for VCF
    pragma_fileformat = generate_custom_unstructured_metainfomation_line("fileformat", "VCFv4.4")
    pragmas_to_add.append(pragma_fileformat)
    print("# start 1")
    print(pragma_fileformat)
    print(pragmas_to_add[-1])
    print("# end 1")
    #TODO:  reference=file, contig, INFO#
    for pragma in GVF_PRAGMAS:
        # file date
        if pragma.startswith("##file-date"):
            date = pragma.split(" ")[1].replace("-", "")
            pragma_filedate = generate_custom_unstructured_metainfomation_line("fileDate", date)
            pragmas_to_add.append(pragma_filedate)
            print("# start 2")
            print(pragma_filedate)
            print(pragmas_to_add[-1])
            print("# end 2")
        # source
        pragma_source = generate_custom_unstructured_metainfomation_line("source", vcf_obj.source)
        pragmas_to_add.append(pragma_source)
    #     # reference #TODO: add this
    #     # contig: recommended, see section 1.4.7 of VCF specification # TODO: add this
    #     # phasing # not required
    #     if pragma.startswith("##gff-version"):
    #         gff_version_number = pragma.split(" ")[1]
    #         pragma_gff_version_number = generate_custom_unstructured_metainfomation_line("gff-version", gff_version_number)
    #         pragmas_to_add.append(pragma_gff_version_number)
    #     elif pragma.startswith("##gvf-version"):
    #         gvf_version_number = pragma.split(" ")[1]
    #         pragma_gvf_version_number = generate_custom_unstructured_metainfomation_line("gvf-version", gvf_version_number)
    #         pragmas_to_add.append(pragma_gvf_version_number)
    #     elif pragma.startswith("##species"):
    #         species_value = pragma.split(" ")[1]
    #         pragma_species_value = generate_custom_unstructured_metainfomation_line("species", species_value)
    #         pragmas_to_add.append(pragma_species_value)
    #     elif pragma.startswith("##genome-build"):
    #         genome_build = pragma.split("genome-build ")[1]
    #         pragma_genome_build = generate_custom_unstructured_metainfomation_line("genome-build", genome_build)
    #         pragmas_to_add.append(pragma_genome_build)
    #     else:
    #         pass
    # for pragma in pragmas_to_add:
    #     if pragma not in unique_pragmas_to_add:
    #         unique_pragmas_to_add.append(pragma)
    # for unique_pragmas in unique_pragmas_to_add:
    #     print("uniq", unique_pragmas)
    # print(all_possible_INFO_lines)
    # print(lines_INFO)
    # #TODO: non-essentials



# step 10
# call step 9
# TODO: finish the below for sample names
def generate_vcf_header(samples):
    vcf_header_fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    for sample in samples:
        vcf_header_fields.append(sample)
    return vcf_header_fields

def generate_vcf_datalines(converted_vcf_lines):
    for line in converted_vcf_lines:
        print(line)



#print("here are the lists:", lines_ALT, lines_INFO, lines_FILTER, lines_FORMAT) # TODO: why is this not populated
#print("here are custom lists", lines_custom_structured, lines_custom_unstructured) # TODO: why is some of this populated.
#TODO: docstrings

# CALLS
VCF_DATA_LINES = {} # DICTIONARY OF LISTS
def main():
    read_in_gvf_file(args)
    read_dgva_info_attributes(dgva_info_attributes_file="dgvaINFOattributes.txt") # needed to generate custom strings
    read_gvf_info_attributes(gvf_info_attributes_file="gvfINFOattributes.txt")

    list_of_vcf_objects = []
    # create a vcf object for every feature line in the GVF (1:1)
    # add the newly created vcf object to the vcf data line it belongs to
    # (1:many; key=chrom_pos; 1 key: many vcf objects)
    for gvf_line_obj in GVF_LINES_OBJ_LIST:
        vcf_object = VcfDataObj(gvf_line_obj)
        list_of_vcf_objects.append(vcf_object)
        if vcf_object.key in VCF_DATA_LINES:
            VCF_DATA_LINES[vcf_object.key].append(vcf_object)
        else:
            vcf_object_dataline_list = []
            vcf_object_dataline_list.append(vcf_object)
            VCF_DATA_LINES[vcf_object.key] = vcf_object_dataline_list
        print("# start - creating vcf object - 1")
        print(VCF_DATA_LINES[vcf_object.key][-1].get_string())
        print("# end - creating vcf object - 1")
        print("number of keys in dict", len(VCF_DATA_LINES.keys()))

        # check the number of objects to see if they are merged
        for key in VCF_DATA_LINES.keys():
            vcf_obj_list = VCF_DATA_LINES[key]
            print("for ", key, " the number of vcf objects is: ", len(vcf_obj_list))
    # # 10c
    # converted_lines = []
    # for vcf_obj in list_of_vcf_objects:
    #     vcf_line = (f"{vcf_obj.chrom}\t"
    #                 f"{vcf_obj.pos}\t"
    #                 f"{vcf_obj.id}\t"
    #                 f"{vcf_obj.ref}\t" #TODO: should this always be empty
    #                 f"{vcf_obj.alt}\t" #TODO: should this always be empty
    #                 f"{vcf_obj.qual}\t" #TODO: should this always be empty
    #                 f"{vcf_obj.filter}\t" #TODO: should this always be empty
    #                 f"{vcf_obj.info}\t"
    #                 f"{vcf_obj.format}\tsampleFORMAThere" #TODO: fill this in
    #                 )
    #     converted_lines.append(vcf_line)
    # print("generating entire VCF header")
    # generate_vcf_metainformation()
    # iterate over meta information line
    # generate standard VCF header
    # sam = ["samA", "samB", "samC"]
    # print(generate_vcf_header(sam))
    # iterate over vcf_dataline list and print vcf data line to output, see 10c
    #generate_vcf_datalines(converted_lines)



if __name__ == "__main__":
    print("Running the GVF to VCF converter")
    # step 1
    parser = argparse.ArgumentParser()
    parser.add_argument("gvf_input", help="GVF input file.")
    parser.add_argument("vcf_output", help="VCF output file.")
    args = parser.parse_args()
    print("The provided input file is: ", args.gvf_input)
    print("The provided output file is: ", args.vcf_output)
    main()
