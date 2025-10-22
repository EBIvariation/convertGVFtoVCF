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
        gvf_attribute_dictionary = get_gvf_attributes(column9_of_gvf)
        vcf_vals = {}
        catching_for_review = []
        # print("dgva_attribute_dict", dgva_attribute_dict)
        mapping_attribute_dict = read_info_attributes(os.path.join(etc_folder, 'attribute_mapper.tsv'))
        # created a rough guide to attributes_for_custom_structured_metainformation in INFOattributes.tsv = this probably should be refined at a later date
        # TODO: edit INFOattributes.tsv i.e. replace unknown placeholders '.' with the actual answer, provide a more informative description
        for attrib_key in gvf_attribute_dictionary:
            # if dgva specific key, create custom INFO tag's meta information line
            if attrib_key in info_attribute_dict:
                field_lines_dictionary["INFO"].append(
                    generate_custom_structured_meta_line(
                        vcf_key="INFO", vcf_key_id=attrib_key,
                        vcf_key_number=info_attribute_dict[attrib_key][1],
                        vcf_key_type=info_attribute_dict[attrib_key][2],
                        vcf_key_description=info_attribute_dict[attrib_key][3],
                        optional_extra_fields=None)
                )
                vcf_vals[attrib_key] = gvf_attribute_dictionary[attrib_key]
            elif attrib_key in mapping_attribute_dict:
                field = mapping_attribute_dict[attrib_key][1]
                key_for_field = mapping_attribute_dict[attrib_key][2]
                field_lines_dictionary[field].append(all_possible_lines_dictionary[field][key_for_field])

            elif attrib_key == "sample_name":
                # sample_names.append(sample_names)
                pass
            # GVF keys (not dgva specific)
            elif attrib_key == "ID":
                pass
            elif attrib_key == "Variant_seq":
                pass
            elif attrib_key == "Reference_seq":
                pass

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
            # elif attrib_key == "Total_reads":
            #     # reserved info key, DP
            #     pass
            # elif attrib_key == "Variant_freq":
            #     # reserve info tag, AF
            #     pass
            # elif attrib_key == "Genotype":
            #     # GT
            #     pass
            else:
                print("catching these attribute keys for review at a later date", attrib_key)
                catching_for_review.append(attrib_key)
        # print("dictionary", gvf_attribute_dictionary)
        # print("vcf_vals", vcf_vals)
        return gvf_attribute_dictionary