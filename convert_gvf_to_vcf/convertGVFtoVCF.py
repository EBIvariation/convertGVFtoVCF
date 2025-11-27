import argparse
import os
from convert_gvf_to_vcf.utils import read_pragma_mapper, read_in_gvf_file
from convert_gvf_to_vcf.vcfline import VcfLine
from convert_gvf_to_vcf.logger import set_up_logging, logger
from convert_gvf_to_vcf.lookup import Lookup
# setting up paths to useful directories
convert_gvf_to_vcf_folder = os.path.dirname(__file__)
etc_folder = os.path.join(convert_gvf_to_vcf_folder, 'etc')

# the functions below relate to the VCF header (Part 1)
def generate_vcf_header_structured_lines(header_type, mapping_attribute_dict):
    """ Generates a dictionary of all possible standard structured lines for INFO/FILTER/FORMAT/ALT.
    :param header_type: type of header file to read i.e. ALT, FILTER, INFO or FORMAT
    :param mapping_attribute_dict: dictionary of all attributes
    :return: dictionary of all possible standard structured lines keys for the header type
    """
    all_possible_lines = {}

    for attribute in mapping_attribute_dict:
        # Formatting the header string for FILTER, INFO or FORMAT and storing in a dictionary
        if mapping_attribute_dict[attribute].get(header_type) is not None and header_type != "ALT":
            header_string = (f'##{header_type}='
                             f'<ID={mapping_attribute_dict[attribute][header_type]["FieldKey"]},'
                             f'Number={mapping_attribute_dict[attribute][header_type]["Number"]},'
                             f'Type={mapping_attribute_dict[attribute][header_type]["Type"]},'
                             f'Description="{mapping_attribute_dict[attribute][header_type]["Description"]}">')
            all_possible_lines[mapping_attribute_dict[attribute][header_type]["FieldKey"]] = header_string
        # Formatting the header string for ALT and storing in a dictionary
        elif mapping_attribute_dict[attribute].get(header_type) is not None and header_type == "ALT":
            if mapping_attribute_dict[attribute][header_type]["FieldKey"] is not None:
                header_string = (f'##{header_type}='
                                 f'<ID={mapping_attribute_dict[attribute][header_type]["FieldKey"]},'
                                 f'Description="{mapping_attribute_dict[attribute][header_type]["Description"]}">')
                all_possible_lines[mapping_attribute_dict[attribute][header_type]["FieldKey"]] = header_string
        else:
            pass
    return all_possible_lines

def generate_vcf_header_unstructured_line(vcf_unstructured_key,
                                          vcf_unstructured_value):
    """ Generates a formatted unstructured metainformation line using a custom key value pair e.g. "##key=value"
    :param vcf_unstructured_key: key for custom unstructured metainformation line
    :param vcf_unstructured_value: value for custom unstructured metainformation line
    :return: custom_unstructured_string
    """
    custom_unstructured_string = f"##{vcf_unstructured_key}={vcf_unstructured_value}"
    return custom_unstructured_string

def generate_vcf_header_metainfo(gvf_pragmas,
                                 gvf_non_essential,
                                 list_of_vcf_objects,
                                 standard_lines_dictionary):
    """ Generates a list of metainformation lines for the VCF header
    :param gvf_pragmas: list of gvf pragmas to convert
    :param gvf_non_essential: list of non-essential gvf pragmas to convert
    :param list_of_vcf_objects: list of vcf objects
    :param standard_lines_dictionary: dictionary of standard lines
    :return: unique_pragmas_to_add, sample_names: a list of pragmas (removed duplicates), list of sample names
    """
    pragmas_to_add = []
    unique_pragmas_to_add = []
    sample_names = []
    unique_alt_lines_to_add = []
    unique_info_lines_to_add = []
    unique_filter_lines_to_add = []
    unique_format_lines_to_add = []
    ####
    # MANDATORY: file format for VCF
    pragmas_to_add.append(generate_vcf_header_unstructured_line("fileformat", "VCFv4.4"))

    #Go through essential pragmas
    #TODO: list of pragmas to add:reference=file, contig, phasing,INFO#
    list_of_pragma = ["##file-date", "##gff-version", "##gvf-version", "##species", "##genome-build"]
    pragma_to_vcf_map = read_pragma_mapper(os.path.join(etc_folder, 'pragma_mapper.tsv'))
    for pragma in gvf_pragmas:
        vcf_header_key, pragma_name, pragma_value = get_pragma_name_and_value(pragma, " ", list_of_pragma, pragma_to_vcf_map)
        pragmas_to_add.append(generate_vcf_header_unstructured_line(vcf_header_key, pragma_value))
        for vcf_obj in list_of_vcf_objects:
            pragmas_to_add.append(generate_vcf_header_unstructured_line("source", vcf_obj.source))
    ####
    ####
    # Go through non-essential pragmas
    list_of_non_essential_pragma = ["#sample", "#Study_accession", "#Study_type", "#Display_name", "#Publication"
                                    "#Study", "#Assembly_name", "#subject"]
    for non_essential_pragma in gvf_non_essential:
        vcf_header_key, pragma_name, pragma_value = get_pragma_name_and_value(non_essential_pragma, ": ", list_of_non_essential_pragma, pragma_to_vcf_map)
        if pragma_name.startswith("#Publication"):
            publication_tokens = get_pragma_tokens(pragma_value, ";", "=")
            pragmas_to_add.append(generate_vcf_header_unstructured_line(publication_tokens[0], publication_tokens[1]))
        elif pragma_name == "#Study":
            study_tokens = get_pragma_tokens(pragma_value, ";", "=")
            pragmas_to_add.append(generate_vcf_header_unstructured_line(study_tokens[0], study_tokens[1]))
        else:
            if vcf_header_key is not None:
                pragmas_to_add.append(generate_vcf_header_unstructured_line(vcf_header_key, pragma_value))
        ####
        # populating sample headers
        if pragma_name.startswith("#sample"):
            list_of_sample_information = pragma_value.split(";")
            for sample_info in list_of_sample_information:
                if sample_info.startswith("sample_name"):
                    sample_names.append(sample_info.split("=")[1])
    # ensure unique sample names and preserve order
    seen_sample_names = set()
    uniq_sample_name = []
    for sample in sample_names:
        if sample not in seen_sample_names:
            seen_sample_names.add(sample)
            uniq_sample_name.append(sample)
    ###
    unique_pragmas_to_add = list(dict.fromkeys(pragma for pragma in pragmas_to_add if pragma not in unique_pragmas_to_add))
    unique_alt_lines_to_add = list(dict.fromkeys(alt_line for alt_line in standard_lines_dictionary["ALT"] if alt_line not in unique_alt_lines_to_add))
    unique_info_lines_to_add = list(dict.fromkeys(info_line for info_line in standard_lines_dictionary["INFO"] if info_line not in unique_info_lines_to_add))
    unique_filter_lines_to_add = list(dict.fromkeys(filter_line for filter_line in standard_lines_dictionary["FILTER"] if filter_line not in unique_filter_lines_to_add))
    unique_format_lines_to_add = list(dict.fromkeys(format_line for format_line in standard_lines_dictionary["FORMAT"] if format_line not in unique_format_lines_to_add))

    return unique_pragmas_to_add, uniq_sample_name, unique_alt_lines_to_add, unique_info_lines_to_add, unique_filter_lines_to_add, unique_format_lines_to_add

# the function below relates to the VCF headerline (Part 2)
def generate_vcf_header_line(samples):
    """ Generates the VCF header line using the nine mandatory headers and the sample names.
    :param samples: list of samples, these will appear in the header line
    :return: vcf_header: a string
    """
    vcf_header_fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    for sample in samples:
        vcf_header_fields.append(sample)
    vcf_header = '\t'.join(vcf_header_fields)
    return vcf_header

# the functions below relate to the GVF header
def parse_pragma(pragma_to_parse, delimiter):
    """ Parses pragma and returns name and value of the pragma.
    :param pragma_to_parse: pragma
    :param delimiter: to split by
    :return: pragma_name, pragma_value: key and value of pragma
    """
    try:
        pragma_tokens = pragma_to_parse.split(delimiter)
        pragma_name = pragma_tokens[0]
        if len(pragma_tokens) >= 2:
            pragma_value = ''.join(map(str, pragma_tokens[1:]))
        # elif len(pragma_tokens) == 1:
            # pragma_value = ''.join(map(str, pragma_tokens[0]))
        else:
            pragma_value = None
            logger.warning(f"WARNING: no value for the following pragma {pragma_to_parse}")
        return pragma_name, pragma_value
    except ValueError:
        logger.error(f"Skipping this, can't be parsed {pragma_to_parse}")

def get_pragma_name_and_value(pragma_to_parse, delimiter, pragma_list, pragma_name_to_vcf_dict):
    """Get pragma name and value and its corresponding VCF header key.
    :param pragma_to_parse: pragma that will be parsed
    :param delimiter: the separator
    :param pragma_list: list of pragmas to search through
    :param pragma_name_to_vcf_dict: dictionary pragma name and its vcf entry
    :return vcf_header_key, pragma_name, pragma_value
    """
    pragma_name, pragma_value = parse_pragma(pragma_to_parse, delimiter)
    if pragma_name in pragma_list:
        vcf_header_key = pragma_name_to_vcf_dict.get(pragma_name)
    else:
        vcf_header_key = None
    return vcf_header_key, pragma_name, pragma_value

def get_pragma_tokens(pragma_value, first_delimiter, second_delimiter):
    """Get pragma tokens for nested pragmas
    :param pragma_value: value to parse
    :param first_delimiter: first separator
    :param second_delimiter: second separtor
    :return pragma_tokens
    """
    initial_list = pragma_value.split(first_delimiter)
    pragma_tokens = []
    for element in initial_list:
        pragma_tokens = element.split(second_delimiter)
    return pragma_tokens

# This is the main conversion logic
def convert_gvf_features_to_vcf_objects(gvf_lines_obj_list, reference_lookup):
    """ Creates VCF objects from GVF feature lines and stores the VCF objects.
    :param gvf_lines_obj_list: list of GVF feature line objects
    :param reference_lookup: an object that stores important dictionaries to be used for reference lookups.
    :return: standard_header_lines, vcf_data_lines, list_of_vcf_objects: header lines for this VCF, datalines for this VCF and a list of VCF objects
    """
    vcf_data_lines = {}  # DICTIONARY OF LISTS, {Chromosome_Pos: [VCF line object]}
    list_of_vcf_objects = []
    # Create data structure to store the header lines for this VCF file (standard meta-information lines)
    standard_header_lines ={
        "ALT": [],
        "INFO": [],
        "FILTER": [],
        "FORMAT": [],
    }
    #TODO: place the all_header_lines_per_type_dict into the reference_lookup.

    # Create data structure to store all possible outcomes for header lines (for fields ALT, INFO, FILTER, FORMAT)
    all_header_lines_per_type_dict = {
        htype: generate_vcf_header_structured_lines(htype, reference_lookup.mapping_attribute_dict) for htype in ["ALT", "INFO", "FILTER", "FORMAT"]
    }

    # Create a vcf object for every feature line in the GVF (1:1)
    for gvf_featureline in gvf_lines_obj_list:
        #NOTE: this is the main Logic of the code
        vcf_object = VcfLine(gvf_featureline,
                             standard_header_lines,
                             all_header_lines_per_type_dict,
                             reference_lookup)
        # Store VCF object in the list
        list_of_vcf_objects.append(vcf_object)

        # vcf_object.key is formatted as follows: Chromosome_Pos
        if vcf_object.key in vcf_data_lines:
            # Add VCF object to the dictionary of lists
            vcf_data_lines[vcf_object.key].append(vcf_object)
        else:
            # Get it into a format where the VCF object can be added to the dictionary of lists
            vcf_data_line_objects_list = [vcf_object]
            vcf_data_lines[vcf_object.key] = vcf_data_line_objects_list
    # Returns the header of the VCF file, the datalines of the VCF file, and the object.
    return standard_header_lines, vcf_data_lines, list_of_vcf_objects

# The functions below relate to the VCF objects
def compare_vcf_objects(list_of_vcf_objects):
    """ Compares VCF objects in the list with the VCF object before it. Returns boolean values.
    :params: list_of_vcf_objects: list of vcf objects
    :return: comparison_results: list of booleans. For future reference, if True, this will determine merging lines; if False, this will determine use of the previous line.
    """
    comparison_results = []
    # For each vcf line object, compare with the previous vcf line object in the list
    for index in range(1, len(list_of_vcf_objects)):
        current_vcf_object = list_of_vcf_objects[index]
        previous_vcf_object = list_of_vcf_objects[index - 1]
        # Determines the VCF line objects as equal based on the CHROM, POS and REF being the same (__eq__ in Vcfline)
        if current_vcf_object == previous_vcf_object:
            comparison_results.append(True) # This will use require merging.
        else:
            comparison_results.append(False) # No merging required. Use previous object.
    return comparison_results

def merge_vcf_objects(previous, current, list_of_sample_names):
    """ Merge VCF objects.
    :params: previous: previous VCF line object
    :params: current: current VCF line object
    :params: list_of_sample_names: sample names
    :return: merged_object
    """
    merged_object = previous.merge(current, list_of_sample_names)
    return merged_object

def keep_vcf_objects(previous, list_of_sample_names):
    """ Keep VCF objects.
    :params: previous VCF line object
    :return: kept_object
    """
    kept_object = previous.keep(list_of_sample_names)
    return kept_object

def determine_merge_or_keep_vcf_objects(list_of_vcf_objects, comparison_results, list_of_sample_names):
    """ Runs through the list of VCF objects and its corresponding comparison result.
    If True, merge parts of the vcf object together. If False, use the previous object
    :params: list_of_vcf_objects: list of vcf line objects
    :return: merge_or_kept_objects: list of vcf line objects that have either been merged or kept as is.
    """
    merge_or_kept_objects = []
    # start at 1 to ensure the first element has a previous object
    for index, compare_result in enumerate(comparison_results, start=1):
        # Merge if the previous and current VCF object are the same (compare_result is True)
        if compare_result:
            merged_object = merge_vcf_objects(list_of_vcf_objects[index - 1], list_of_vcf_objects[index], list_of_sample_names)
            merge_or_kept_objects.append(merged_object)
        # Keep previous if previous and current VCF object are different (compare_result is False)
        else:
            # Prevents duplicates: Ensure previous line is not the same as the line before it.
            if list_of_vcf_objects[index - 2] != list_of_vcf_objects[index - 1 ]:
                # keep the previous VCF line object
                kept_object = keep_vcf_objects(list_of_vcf_objects[index - 1], list_of_sample_names)
                merge_or_kept_objects.append(kept_object)
    return merge_or_kept_objects

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("gvf_input", help="GVF input file.")
    parser.add_argument("vcf_output", help="VCF output file.")
    parser.add_argument("-a", "--assembly", help="FASTA assembly file")
    parser.add_argument("--log", help="Path to log file")
    args = parser.parse_args()

    # Set up logging functionality
    if args.log:
        log_path = set_up_logging(args.log)
    else:
        log_path = set_up_logging()

    # Log the inputs and outputs.
    logger.info("Running the GVF to VCF converter")
    logger.info(f"The provided input file is: {args.gvf_input}")
    logger.info(f"The provided output file is: {args.vcf_output}")
    if args.assembly:
        logger.info(f"The provided assembly file is: {args.assembly}")
    assembly_file = os.path.abspath(args.assembly)
    assert os.path.isfile(assembly_file), "Assembly file does not exist"
    logger.info(f"The log file is {log_path}")

    # Read input file and separate out its components
    logger.info(f"Reading in the following GVF input: {args.gvf_input}")
    gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(args.gvf_input)

    # Creating lookup object to store important dictionaries and log what has been stored.
    reference_lookup = Lookup(assembly_file)
    logger.info("Creating the reference lookup object.")
    logger.info("Storing the attributes file: attribute_mapper.yaml")
    logger.info("Storing the symbolic allele dictionary.")
    logger.info(f"Storing the assembly file: {assembly_file}")
    logger.info("Storing the IUPAC ambiguity dictionary.")

    # Convert each feature line in the GVF file to a VCF object (stores all the data for a line in the VCF file).
    # NOTE: Main Logic lives here.
    (
        header_lines,
        vcf_data_lines, #TODO: check if this can be removed
        list_of_vcf_objects
    ) = convert_gvf_features_to_vcf_objects(gvf_lines_obj_list, reference_lookup)

    logger.info(f"Writing to the following VCF output: {args.vcf_output}")
    logger.info("Generating the VCF header and the meta-information lines")
    with open(args.vcf_output, "w") as vcf_output:
        # Preparation work:
        # Store the VCF metainformation and ensure preservation of important GVF data.
        # This information will be useful when creating the VCF header.
        # TODO: refactor function generate_vcf_metainfo
        (
            unique_pragmas_to_add,
            samples,
            unique_alt_lines_to_add,
            unique_info_lines_to_add,
            unique_filter_lines_to_add,
            unique_format_lines_to_add
        ) = generate_vcf_header_metainfo(gvf_pragmas,
                                         gvf_non_essential,
                                         list_of_vcf_objects,
                                         header_lines)
        logger.info(f"Total number of samples in this VCF: {len(samples)}")

        # Part 1 of VCF file: Write the VCF header. This will include perserved data from the GVF file.
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

        # Part 2 of VCF file: Write the VCF header line. This is the nine mandatory fields with its sample names.
        header_fields = generate_vcf_header_line(samples)
        vcf_output.write(f"{header_fields}\n")

        # Part 3 of VCF file: Write the VCF data lines. This will contain info about the position in the genome,
        # its variants and genotype information per sample.
        logger.info("Generating the VCF datalines")
        # Each GVF feature has been converted to a VCF object so begin comparing and merging the VCF objects.
        comparison_flags = compare_vcf_objects(list_of_vcf_objects) # Identifies which VCF objects to merge
        merge_or_kept_vcf_objects = determine_merge_or_keep_vcf_objects(list_of_vcf_objects, comparison_flags, samples)
        # Write the VCF objects as data lines in the VCF file.
        for vcf_line_object in merge_or_kept_vcf_objects:
            vcf_output.write(str(vcf_line_object) + "\n")
            # vcf_output.write("\t".join(str(val) for val in line) + "\n")
    vcf_output.close()
    logger.info("GVF to VCF conversion complete")


if __name__ == "__main__":
    main()
