import argparse
import os
from convert_gvf_to_vcf.utils import read_pragma_mapper, \
    read_in_gvf_file, \
    read_yaml, generate_symbolic_allele_dict
from convert_gvf_to_vcf.vcfline import VcfLine
from convert_gvf_to_vcf.logger import set_up_logging, logger

# setting up paths to useful directories
convert_gvf_to_vcf_folder = os.path.dirname(__file__)
etc_folder = os.path.join(convert_gvf_to_vcf_folder, 'etc')

def generate_vcf_header_structured_lines(header_type, mapping_attribute_dict):
    """ Generates a dictionary of all possible standard structured lines for INFO/FILTER/FORMAT/ALT
    :param header_type: type of header file to read i.e. ALT, FILTER, INFO or FORMAT
    :param mapping_attribute_dict: dictionary of all attributes
    :return: dictionary of all possible standard structured lines keys for the header type
    """
    all_possible_lines = {}

    for attribute in mapping_attribute_dict:
        if mapping_attribute_dict[attribute].get(header_type) is not None and header_type != "ALT":
            header_string = (f'##{header_type}='
                             f'<ID={mapping_attribute_dict[attribute][header_type]["FieldKey"]},'
                             f'Number={mapping_attribute_dict[attribute][header_type]["Number"]},'
                             f'Type={mapping_attribute_dict[attribute][header_type]["Type"]},'
                             f'Description="{mapping_attribute_dict[attribute][header_type]["Description"]}">')
            all_possible_lines[mapping_attribute_dict[attribute][header_type]["FieldKey"]] = header_string
        elif mapping_attribute_dict[attribute].get(header_type) is not None and header_type == "ALT":
            if mapping_attribute_dict[attribute][header_type]["FieldKey"] is not None:
                header_string = (f'##{header_type}='
                                 f'<ID={mapping_attribute_dict[attribute][header_type]["FieldKey"]},'
                                 f'Description="{mapping_attribute_dict[attribute][header_type]["Description"]}">')
                all_possible_lines[mapping_attribute_dict[attribute][header_type]["FieldKey"]] = header_string
        else:
            pass
    return all_possible_lines

def generate_custom_unstructured_meta_line(vcf_unstructured_key,
                                           vcf_unstructured_value):
    """ Generates a formatted unstructured metainformation line using a custom key value pair.
    :param vcf_unstructured_key: key for custom unstructured metainformation line
    :param vcf_unstructured_value: value for custom unstructured metainformation line
    :return: custom_unstructured_string
    """
    custom_unstructured_string = f"##{vcf_unstructured_key}={vcf_unstructured_value}"
    return custom_unstructured_string


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
            logger.warning("WARNING: no value for the following pragma %s", pragma_to_parse)
        return pragma_name, pragma_value
    except ValueError:
        logger.error("Skipping this, can't be parsed %s", pragma_to_parse)

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

def generate_vcf_metainfo(gvf_pragmas, gvf_non_essential, list_of_vcf_objects,
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
    # MANDATORY: file format for VCF
    pragmas_to_add.append(generate_custom_unstructured_meta_line("fileformat", "VCFv4.4"))
    #Go through essential pragmas
    #TODO: list of pragmas to add:reference=file, contig, phasing,INFO#
    list_of_pragma = ["##file-date", "##gff-version", "##gvf-version", "##species", "##genome-build"]
    pragma_to_vcf_map = read_pragma_mapper(os.path.join(etc_folder, 'pragma_mapper.tsv'))
    for pragma in gvf_pragmas:
        vcf_header_key, pragma_name, pragma_value = get_pragma_name_and_value(pragma, " ", list_of_pragma, pragma_to_vcf_map)
        pragmas_to_add.append(generate_custom_unstructured_meta_line(vcf_header_key, pragma_value))
        for vcf_obj in list_of_vcf_objects:
            pragmas_to_add.append(generate_custom_unstructured_meta_line("source", vcf_obj.source))
    # Go through non-essential pragmas
    list_of_non_essential_pragma = ["#sample", "#Study_accession", "#Study_type", "#Display_name", "#Publication"
                                    "#Study", "#Assembly_name", "#subject"]
    for non_essential_pragma in gvf_non_essential:
        vcf_header_key, pragma_name, pragma_value = get_pragma_name_and_value(non_essential_pragma, ": ", list_of_non_essential_pragma, pragma_to_vcf_map)
        if pragma_name.startswith("#Publication"):
            publication_tokens = get_pragma_tokens(pragma_value, ";", "=")
            pragmas_to_add.append(generate_custom_unstructured_meta_line(publication_tokens[0], publication_tokens[1]))
        elif pragma_name == "#Study":
            study_tokens = get_pragma_tokens(pragma_value, ";", "=")
            pragmas_to_add.append(generate_custom_unstructured_meta_line(study_tokens[0], study_tokens[1]))
        else:
            if vcf_header_key is not None:
                pragmas_to_add.append(generate_custom_unstructured_meta_line(vcf_header_key, pragma_value))
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

    unique_pragmas_to_add = list(dict.fromkeys(pragma for pragma in pragmas_to_add if pragma not in unique_pragmas_to_add))
    unique_alt_lines_to_add = list(dict.fromkeys(alt_line for alt_line in standard_lines_dictionary["ALT"] if alt_line not in unique_alt_lines_to_add))
    unique_info_lines_to_add = list(dict.fromkeys(info_line for info_line in standard_lines_dictionary["INFO"] if info_line not in unique_info_lines_to_add))
    unique_filter_lines_to_add = list(dict.fromkeys(filter_line for filter_line in standard_lines_dictionary["FILTER"] if filter_line not in unique_filter_lines_to_add))
    unique_format_lines_to_add = list(dict.fromkeys(format_line for format_line in standard_lines_dictionary["FORMAT"] if format_line not in unique_format_lines_to_add))

    return unique_pragmas_to_add, uniq_sample_name, unique_alt_lines_to_add, unique_info_lines_to_add, unique_filter_lines_to_add, unique_format_lines_to_add

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
                                assembly_file, mapping_attribute_dict,
                                symbolic_allele_dictionary
                                ):
    """ Creates VCF objects from GVF feature lines and stores the VCF objects.
    :param gvf_lines_obj_list: list of GVF feature line objects
    :param assembly_file: FASTA file to assembly
    :param mapping_attribute_dict: dictionary of attributes
    :param symbolic_allele_dictionary: symbolic_allele_dictionary
    :return: standard_header_lines, vcf_data_lines, list_of_vcf_objects: header lines for this VCF, datalines for this VCF and a list of VCF objects
    """
    vcf_data_lines = {}  # DICTIONARY OF LISTS
    list_of_vcf_objects = []
    # standard meta-information lines for this VCF file
    standard_header_lines ={
        "ALT": [],
        "INFO": [],
        "FILTER": [],
        "FORMAT": [],
    }
    all_header_lines_per_type_dict = {
        htype: generate_vcf_header_structured_lines(htype, mapping_attribute_dict) for htype in ["ALT", "INFO", "FILTER", "FORMAT"]
    }

    # create a vcf object for every feature line in the GVF (1:1)
    # add the newly created vcf object to the vcf data line it belongs to
    # (1:many; key=chrom_pos; 1 key: many vcf objects)
    for gvf_featureline in gvf_lines_obj_list:
        vcf_object = VcfLine(gvf_featureline,
                             mapping_attribute_dict,
                             symbolic_allele_dictionary,
                             assembly_file,
                             standard_header_lines,
                             all_header_lines_per_type_dict)


        list_of_vcf_objects.append(vcf_object)
        if vcf_object.key in vcf_data_lines:
            vcf_data_lines[vcf_object.key].append(vcf_object)
        else:
            vcf_data_line_objects_list = [vcf_object]
            vcf_data_lines[vcf_object.key] = vcf_data_line_objects_list
        # check the number of objects to see if they are merged
        # for key in vcf_data_lines.keys():
        #     vcf_obj_list = vcf_data_lines[key]
        #     print("for", key, " the number of vcf objects is: ", len(vcf_obj_list))
    return standard_header_lines, vcf_data_lines, list_of_vcf_objects

def format_sample_values(sample_name_dict_format_kv, list_of_sample_names):
    """ Creates a partial vcf data line of sample format values.
    :param sample_name_dict_format_kv: dictionary of sample names => sample format value
    :param list_of_sample_names: list of sample names
    :return: sample_format_values_string: formatted string
    """
    sample_format_value_tokens = []
    for sample in list_of_sample_names:
        if sample in sample_name_dict_format_kv:
            # only one format value
            format_value = sample_name_dict_format_kv[sample]
            if len(format_value) == 1:
                for onesample in list_of_sample_names:
                    if onesample != sample:
                        sample_format_value_tokens.append(".")
                    else:
                        sample_format_value_tokens.append(':'.join(format_value.values()))
            elif len(format_value) > 1:  # more than one format value
                # if the sample is found, go through all other sample names, and add the key with value '.'
                keys_to_populate = format_value.keys()
                for sam in list_of_sample_names:
                    if sam != sample:
                        sample_name_dict_format_kv[sam] = {}
                        for key_to_populate in keys_to_populate:
                            sample_name_dict_format_kv[sam][key_to_populate] = "."
                    multi_format_value = sample_name_dict_format_kv[sam]
                    sample_format_value_tokens.append(':'.join(multi_format_value.values()))
            else: # in the event of an empty format column
                format_value = "." # set to missing value
                sample_format_value_tokens.append(format_value)
        else:
            all_missing = all(sample not in sample_name_dict_format_kv for sample in list_of_sample_names)
    # deal with empty format tag
    if all_missing:
        format_value = "."
        for s in list_of_sample_names:
            sample_format_value_tokens.append(format_value)

    sample_format_values_string = '\t'.join(sample_format_value_tokens)
    return sample_format_values_string

def format_vcf_datalines(list_of_vcf_objects, list_of_sample_names):
    """ Iterates through a list of VCF objects and sample names and formats them as a VCF dataline.
    :param list_of_vcf_objects: list of vcf objects
    :param list_of_sample_names: list of sample names
    :return: formatted_vcf_datalines: list of formatted vcf datalines
    """
    formatted_vcf_datalines = []
    for vcf_obj in list_of_vcf_objects:
        sample_name_dict_format_kv = vcf_obj.format_dict
        sample_format_values_string = format_sample_values(sample_name_dict_format_kv, list_of_sample_names)
        vcf_info_string = ";".join([inf for inf in vcf_obj.info if inf is not None])
        vcf_line = (f"{vcf_obj.chrom}\t"
                        f"{vcf_obj.pos}\t"
                        f"{vcf_obj.id}\t"
                        f"{vcf_obj.ref}\t"
                        f"{vcf_obj.alt}\t"
                        f"{vcf_obj.qual}\t"
                        f"{vcf_obj.filter}\t"
                        #f"{vcf_obj.info}\t"
                        f"{vcf_info_string}\t"
                        f"{vcf_obj.format}\t"
                        f"{sample_format_values_string}"
                        )
        formatted_vcf_datalines.append(vcf_line)
    return formatted_vcf_datalines

def get_bigger_dictionary(dict1, dict2):
    """Determines the biggest of two dictionaries
    :param: dictionary1
    :param: dictinary2
    :return: smallest, largest
    """
    if len(dict1) > len(dict2):
        biggest_dict = dict1
        smallest_dict = dict2
    elif len(dict1) < len(dict2):
        biggest_dict = dict2
        smallest_dict = dict1
    else:
        biggest_dict = dict1
        smallest_dict = dict2
    return smallest_dict, biggest_dict

def merge_and_add(previous_element, current_element, delimiter):
    """ If same, use current element. If different, merge with delimiter.
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

def compare_and_merge_lines(list_of_formatted_vcf_datalines, headerline):
    merged_lines = []
    for previous, current in zip(list_of_formatted_vcf_datalines, list_of_formatted_vcf_datalines[1:]):
        # print(f"previous line:\n{previous}\ncurrent line:\n{current}\n")
        previous_tokens = previous.split("\t")
        current_tokens = current.split("\t")
        header_fields = headerline.split("\t")

        previous_data = dict(zip(header_fields, previous_tokens))
        current_data = dict(zip(header_fields, current_tokens))
        merged_data = {}
        if (
                previous_data["#CHROM"] == current_data["#CHROM"]
                and previous_data["POS"] == current_data["POS"]
                and previous_data["REF"] == current_data["REF"]
        ):
            print("True - merge")
            merged_data["#CHROM"] = current_data["#CHROM"]
            merged_data["POS"] = current_data["POS"]
            merged_data["ID"] = merge_and_add(previous_data["ID"], current_data["ID"], ";")
            merged_data["REF"] = current_data["REF"]
            merged_data["ALT"] = merge_and_add(previous_data["ALT"], current_data["ALT"], ",")
            merged_data["QUAL"] = current_data["QUAL"]
            merged_data["FILTER"] = merge_and_add(previous_data["FILTER"], current_data["FILTER"], ";")
            # INFO
            previous_info_tokens = previous_data["INFO"].split(";")
            current_info_tokens = current_data["INFO"].split(";")
            previous_info_dict = {}
            current_info_dict = {}
            for p in previous_info_tokens:
                p_key, p_value = p.split("=")
                previous_info_dict[p_key] = p_value
            for c in current_info_tokens:
                c_key, c_value = c.split("=")
                current_info_dict[c_key] = c_value

            smallest_info_dict, biggest_info_dict = get_bigger_dictionary(previous_info_dict, current_info_dict)

            merged_dictionary = {}

            for key, value in biggest_info_dict.items():
                merged_dictionary.setdefault(key, []).append(value)
                if smallest_info_dict.get(key) != value:
                    merged_dictionary.setdefault(key, []).append(smallest_info_dict.get(key))

            parts_of_info_string = []
            for merge_key, merge_value in merged_dictionary.items():
                cleaned_merge_value = [value for value in merge_value if value is not None]
                conjoined_merge_value = ",".join(cleaned_merge_value)
                part_string = f"{merge_key}={conjoined_merge_value}"
                parts_of_info_string.append(part_string)
            merged_info_string = ';'.join(parts_of_info_string)
            merged_data["INFO"] = merged_info_string
            # FORMAT
            previous_format_key_tokens = previous_data["FORMAT"].split(":")
            current_format_key_tokens = current_data["FORMAT"].split(":")
            merged_format_key_tokens = []
            for format_key in previous_format_key_tokens + current_format_key_tokens:
                if format_key not in merged_format_key_tokens:
                    merged_format_key_tokens.append(format_key)

            merged_format_key_string = ':'.join(merged_format_key_tokens)
            merged_data["FORMAT"] = merged_format_key_string
            # sample values.
            sample_names = header_fields[9:]
            merged_sample_format = {}
            for sample_name in sample_names:
                previous_sample_format_value = dict(zip(previous_data["FORMAT"].split(":"),previous_data[sample_name].split(":")))
                current_sample_format_value = dict(zip(current_data["FORMAT"].split(":"),current_data[sample_name].split(":")))
                smallest_format_dict, biggest_format_dict = get_bigger_dictionary(previous_sample_format_value, current_sample_format_value)
                for k in biggest_format_dict:
                    biggest_value = biggest_format_dict.get(k)
                    smallest_value = smallest_format_dict.get(k)

                    if biggest_value is None or biggest_value == ".":
                        biggest_value = ""
                    if smallest_value is None or smallest_value == ".":
                        smallest_value = ""
                    element = merge_and_add(biggest_value, smallest_value, "")
                    if element == "":
                        element = "."

                    # merged_sample_format[sample_name] = {k: element}
                    if sample_name not in merged_sample_format:
                        merged_sample_format[sample_name] = {}
                        merged_sample_format[sample_name].setdefault(k, []).append(element)
                    else:
                        merged_sample_format[sample_name].setdefault(k, []).append(element)


            values = []
            for sample_name in sample_names:
                for key in merged_sample_format[sample_name]:
                    values.append(merged_sample_format[sample_name][key])
                flat_values = [v2 for v1 in values for v2 in v1 ]
                sample_format_string =':'.join(flat_values)

                merged_data[sample_name] = sample_format_string
            merged_lines.append(merged_data)
            print("---")
        else:
            print("False - keep previous")
            merged_data["#CHROM"] = previous_data["#CHROM"]
            merged_data["POS"] = previous_data["POS"]
            merged_data["ID"] = previous_data["ID"]
            merged_data["REF"] = previous_data["REF"]
            merged_data["ALT"] = previous_data["ALT"]
            merged_data["QUAL"] = previous_data["QUAL"]
            merged_data["FILTER"] = previous_data["FILTER"]
            merged_data["INFO"] = previous_data["INFO"]
            merged_data["FORMAT"] = previous_data["FORMAT"]
            sample_names = header_fields[9:]
            for sample in sample_names:
                merged_data[sample] = previous_data[sample]

            merged_lines.append(merged_data)
            print("---")
    return merged_lines


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("gvf_input", help="GVF input file.")
    parser.add_argument("vcf_output", help="VCF output file.")
    parser.add_argument("-a", "--assembly", help="FASTA assembly file")
    parser.add_argument("--log", help="Path to log file")
    args = parser.parse_args()

    if args.log:
        log_path = set_up_logging(args.log)
    else:
        log_path = set_up_logging()

    logger.info("Running the GVF to VCF converter")
    logger.info(f"The provided input file is: {args.gvf_input}")
    logger.info(f"The provided output file is: {args.vcf_output}")
    if args.assembly:
        logger.info(f"The provided assembly file is: {args.assembly}")
    assembly_file = os.path.abspath(args.assembly)
    assert os.path.isfile(assembly_file), "Assembly file does not exist"
    logger.info(f"The log file is {log_path}")

    # custom meta-information lines for this VCF file
    logger.info(f"Reading in the following GVF input: {args.gvf_input}")
    gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(args.gvf_input)
    # store attributes and symbolic alleles
    mapping_attribute_dict = read_yaml(os.path.join(etc_folder, "attribute_mapper.yaml"))
    logger.info("Reading in the attributes file: " + "attribute_mapper.yaml")
    symbolic_allele_dictionary = generate_symbolic_allele_dict(mapping_attribute_dict)

    (
        header_lines,
        vcf_data_lines,
        list_of_vcf_objects
    ) = gvf_features_to_vcf_objects(gvf_lines_obj_list,
                                    assembly_file,
                                    mapping_attribute_dict,
                                    symbolic_allele_dictionary
                                    )

    logger.info(f"Writing to the following VCF output: {args.vcf_output}")
    logger.info("Generating the VCF header and the meta-information lines")
    with open(args.vcf_output, "w") as vcf_output:
        (
            unique_pragmas_to_add,
            samples,
            unique_alt_lines_to_add,
            unique_info_lines_to_add,
            unique_filter_lines_to_add,
            unique_format_lines_to_add
        ) = generate_vcf_metainfo(gvf_pragmas, gvf_non_essential, list_of_vcf_objects, header_lines)
        logger.info(f"Total number of samples in this VCF: {len(samples)}")
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
        logger.info("Generating the VCF datalines")
        formatted_vcf_datalines = format_vcf_datalines(list_of_vcf_objects, samples)
        merged_lines = compare_and_merge_lines(formatted_vcf_datalines, header_fields)
        for line in merged_lines:
            vcf_output.write("\t".join(str(val) for val in line.values()) + "\n")
    vcf_output.close()
    logger.info("GVF to VCF conversion complete")


if __name__ == "__main__":
    main()
