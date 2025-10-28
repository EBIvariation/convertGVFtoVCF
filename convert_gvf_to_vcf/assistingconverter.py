# this is an assistant converter to help convert gvf attributes
import os
from convert_gvf_to_vcf.utils import read_info_attributes
from convert_gvf_to_vcf.helpers import generate_custom_structured_meta_line
# setting up paths to useful directories
convert_gvf_to_vcf_folder = os.path.dirname(__file__)
etc_folder = os.path.join(convert_gvf_to_vcf_folder, 'etc')

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

class Assistingconverter:
    @staticmethod
    def convert_gvf_attributes_to_vcf_values(column9_of_gvf,
                                             info_attribute_dict,
                                             field_lines_dictionary,
                                             all_possible_lines_dictionary):
        # this converts GVF attributes to a dictionary that will make VCF values
        # this also populates ALT INFO FILTER FORMAT with the correct VCF values.
        gvf_attribute_dictionary = get_gvf_attributes(column9_of_gvf)
        vcf_info_values = {} # key is info field value; value is value
        catching_for_review = []
        mapping_attribute_dict = read_info_attributes(os.path.join(etc_folder, 'attribute_mapper.tsv'))
        print("dgva_attribute_dict", gvf_attribute_dictionary) #GVF line in dgva
        print("info_attribute_dict", info_attribute_dict) # INFOattributes.tsv
        print("field_lines_dictionary", field_lines_dictionary) # dictionaries for ALT INFO FILTER FORMAT
        print("mapping_attribute_dict", mapping_attribute_dict) # attribute mapper
        # created a rough guide to attributes_for_custom_structured_metainformation in INFOattributes.tsv = this probably should be refined at a later date
        # TODO: edit INFOattributes.tsv i.e. replace unknown placeholders '.' with the actual answer, provide a more informative description
        # for each attribute in the DGVa GVF file,
        # if the attributes is in the INFOattributes.tsv add to INFO dictionary
        # else if attirbute is in attribute_mapper.tsv add to its specific dictionary (INFO/FORMAT)
        # else if attribute is a specific case, deal with appropriately.
        for attrib_key in gvf_attribute_dictionary:
            # if dgva specific key, create custom INFO tag's meta information line
            if attrib_key in info_attribute_dict:
                # create and store header line
                field_lines_dictionary["INFO"].append(
                    generate_custom_structured_meta_line(
                        vcf_key="INFO", vcf_key_id=attrib_key,
                        vcf_key_number=info_attribute_dict[attrib_key][1],
                        vcf_key_type=info_attribute_dict[attrib_key][2],
                        vcf_key_description=info_attribute_dict[attrib_key][3],
                        optional_extra_fields=None)
                )
                # store INFO values
                vcf_info_values[attrib_key] = gvf_attribute_dictionary[attrib_key]
            elif attrib_key in mapping_attribute_dict:
                # create and store header line
                field_lines_dictionary[mapping_attribute_dict[attrib_key][1]].append(all_possible_lines_dictionary[mapping_attribute_dict[attrib_key][1]][mapping_attribute_dict[attrib_key][2]])
                # store INFO values
                vcf_info_values[mapping_attribute_dict[attrib_key][2]] = gvf_attribute_dictionary[attrib_key]
            # GVF keys (not dgva specific)
            elif attrib_key == "Dbxref":
                # custom info tag + pase and add to id?
                pass
            elif attrib_key == "Variant_reads":
                # reserved info/format key, AD/AC
                pass
            elif attrib_key == "Zygosity":
                # format and GT tag
                pass
            elif attrib_key == "Phased":
                # GT or FORMAT PS
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

        info_string = ''.join(f'{key}={value};' for key, value in vcf_info_values.items()).rstrip(';')
        # mapping attribute dictionary
        return gvf_attribute_dictionary, info_string