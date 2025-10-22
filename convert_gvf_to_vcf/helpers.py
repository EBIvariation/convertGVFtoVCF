# helpers.py is to prevent circular imports
def generate_custom_structured_meta_line(vcf_key, vcf_key_id, vcf_key_number, vcf_key_type, vcf_key_description,
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
    custom_structured_string = (f'##{vcf_key}=<'
                                f'ID={vcf_key_id},'
                                f'Number={vcf_key_number},'
                                f'Type={vcf_key_type},'
                                f'Description="{vcf_key_description}"'
                                f'{vcf_key_extra_keys}>')
    return custom_structured_string