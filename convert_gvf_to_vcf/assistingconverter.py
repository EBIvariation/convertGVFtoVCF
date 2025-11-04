"""
This is an assistant converter to help convert gvf attributes
"""
import os
from convert_gvf_to_vcf.utils import read_yaml

# setting up paths to useful directories
convert_gvf_to_vcf_folder = os.path.dirname(__file__)
etc_folder = os.path.join(convert_gvf_to_vcf_folder, 'etc')

def generate_custom_structured_meta_line(vcf_key, vcf_key_id, vcf_key_number,
                                         vcf_key_type, vcf_key_description,
                                         optional_extra_fields=None):
    """ Generates a custom structured meta-information line for INFO/FILTER/FORMAT/ALT
    :param vcf_key: required field INFO, FILTER, FORMAT, ALT
    :param vcf_key_id: required field for structured lines ID
    :param vcf_key_number: Number of values included or special character: A or R or G or .
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
    custom_structured_string = (f'##{vcf_key}=<'
                                f'ID={vcf_key_id},'
                                f'Number={vcf_key_number},'
                                f'Type={vcf_key_type},'
                                f'Description="{vcf_key_description}"'
                                f'{vcf_key_extra_keys}>')
    return custom_structured_string

def get_gvf_attributes(column9_of_gvf):
    """Get a dictionary of GVF attributes
    :param column9_of_gvf:  column - the final column of the GVF file
    :return: gvf_attribute_dictionary: a dictionary of attribute keys and their values
    """
    gvf_attribute_dictionary = {}  # attribute key => value
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


def convert_gvf_attributes_to_vcf_values(column9_of_gvf,
                                         field_lines_dictionary,
                                         all_possible_lines_dictionary):
    """Converts GVF attributes to a dictionary that will store VCF values.
    Populates ALT INFO FILTER FORMAT with the correct VCF values.
    :param column9_of_gvf: attributes column of gvf file
    :param field_lines_dictionary: dictionaries for ALT INFO FILTER and FORMAT
    :param all_possible_lines_dictionary: all possible VCF header lines
    :return gvf_attribute_dictionary, info_string: dict of GVF attributes and formatted info string.
    """
    # this converts GVF attributes to a dictionary that will make VCF values
    # this also populates ALT INFO FILTER FORMAT with the correct VCF values.
    gvf_attribute_dictionary = get_gvf_attributes(column9_of_gvf)
    vcf_info_values = {} # key is info field value; value is value
    vcf_format_values = {} # key is format field value; value is value
    catching_for_review = []
    mapping_attribute_dict = read_yaml(os.path.join(etc_folder, 'attribute_mapper.yaml'))

    for attrib_key, attrib_value in gvf_attribute_dictionary.items():

        if attrib_key in mapping_attribute_dict:
            field_values = mapping_attribute_dict[attrib_key]
            # INFO: create and store header line then store value
            field = "INFO"
            if field in field_values:
                field_key = field_values[field]["FieldKey"]
                header = generate_custom_structured_meta_line(
                    vcf_key=field,
                    vcf_key_id=field_values[field]["FieldKey"],
                    vcf_key_number=field_values[field]["Number"],
                    vcf_key_type=field_values[field]["Type"],
                    vcf_key_description=field_values[field]["Description"],
                    optional_extra_fields=None)
                field_lines_dictionary[field].append(header)
                vcf_info_values[field_key] = gvf_attribute_dictionary[attrib_key]
            # FORMAT: create and store header line then store value
            field = "FORMAT"
            if field in field_values:
                field_key = field_values[field]["FieldKey"]
                field_lines_dictionary[field].append(all_possible_lines_dictionary[field][field_key])
                format_key = field_key
                format_value = gvf_attribute_dictionary[attrib_key]
                sample_name = gvf_attribute_dictionary.get("sample_name")
                if sample_name in vcf_format_values:
                    vcf_format_values[sample_name].update({format_key: format_value})
                else:
                    vcf_format_values[sample_name] = {format_key: format_value}
        else:
            print("catching attribute keys for review at a later date", attrib_key, attrib_value)
            catching_for_review.append(attrib_key)
    info_string = ''.join(f'{key}={value};' for key, value in vcf_info_values.items()).rstrip(';')

    return gvf_attribute_dictionary, info_string, vcf_format_values
