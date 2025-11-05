# this file contains readers
import os
import yaml


from convert_gvf_to_vcf.gvffeature import GvfFeatureline

# setting up paths to useful directories
convert_gvf_to_vcf_folder = os.path.dirname(__file__)
etc_folder = os.path.join(convert_gvf_to_vcf_folder, 'etc')


def read_yaml(yaml_file):
    """Reads a yaml file (of attributes) and returns dictionary
    :param yaml_file: file of attributes
    :return: : A dictionary of attributes.
    """
    with open(yaml_file) as y_file:
        attributes = yaml.safe_load(y_file)
    return attributes

def read_pragma_mapper(pragma_mapper_file):
    """ Reads in the pragma mapper file and stores as dictionary.
    :param pragma_mapper_file: A file containing the pragma and their VCF equivalent
    :return: pragma_to_vcf_header: dictionary of pragma => vcf header
    """
    pragma_to_vcf_header = {}
    with open(pragma_mapper_file) as pragma_file:
        next(pragma_file)
        for line in pragma_file:
            pragma_tokens = line.rstrip().split("\t")
            pragma = pragma_tokens[0]
            vcf_header = pragma_tokens[1]
            pragma_to_vcf_header[pragma] = vcf_header
    return pragma_to_vcf_header

def generate_symbolic_allele_dict(mapping_dictionary):
    """Reads in mapping dictionary and returns a symbolic allele dictionary.
    :param mapping_dictionary: mapping dictionary
    :return symbolic_allele_dict: stores information for a symbolic allele
    """
    symbolic_allele_dict = {}
    for attribute in mapping_dictionary:
        header_type= "ALT"
        if mapping_dictionary[attribute].get(header_type) is not None:
            if mapping_dictionary[attribute].get(header_type).get("FieldKey") is not None:
                name = attribute
                sequence_ontology_id = mapping_dictionary[attribute].get(header_type).get("SOID")
                symb_allele = mapping_dictionary[attribute].get(header_type).get("FieldKey")
                description = mapping_dictionary[attribute].get(header_type).get("Description")

                symbolic_allele_dict.setdefault(name, []).append(sequence_ontology_id)
                symbolic_allele_dict.setdefault(name, []).append(symb_allele)
                symbolic_allele_dict.setdefault(name, []).append(description)
    return symbolic_allele_dict


def read_in_gvf_file(gvf_input):
    """ Reads in the user provided GVF file.
    :param gvf_input: arguments.gvf_input
    :return: gvf_pragmas, gvf_non_essential, gvf_lines_obj_list
    """
    gvf_pragmas = []  # list of pragma lines starting with: ##
    gvf_non_essential = []  # list of non-essential lines starting with: #
    features = []
    gvf_lines_obj_list = []  # list of objects when reading in gvf files, one object represents a gvf line

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