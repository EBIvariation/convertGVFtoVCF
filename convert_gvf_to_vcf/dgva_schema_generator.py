import argparse

import pandas as pd
import json


#TODO: ONCE THE JSON SCHEMA HAS BEEN AGREED FINALISED, REMOVE THIS SCRIPT AND THE SPREADSHEET FROM THE REPO.
# THIS LEAVES THE JSON SCHEMA AS THE ONLY SOURCE OF TRUTH.

class SchemaCreator:
    """The responsibility of this class is to create a JSON schema from a tsv file"""
    def __init__(self, input_file_path):
        # For converting metadata to json string.
        self.input_path = input_file_path
        self.data_frame = self.clean_data(self.read_file(self.input_path))
        self.data_string = self.generate_json_string(self.data_frame)
        # JSON SCHEMA VARIABLES
        self.schema_version = "http://json-schema.org/draft-07/schema#"
        self.version = "1.0.0"
        self.id = "dgva_additional_metadata"
        self.author = "EVA"
        self.title = "Metadata for Database of Genomic Variants archive (DGVa) that are not already captured in EVA metadata schema."
        self.description = "This contains metadata originally found in DGVa and which could not be imported into the European Variation Archive (EVA)."


    def read_file(self, file_path):
        """Reads in table/spreadsheet of Metadata
        :param file_path: path to tsv file of metadata
        :returns: pandas dataframe
        """
        df = pd.read_csv(file_path, sep="\t")
        return df

    def map_type_oracle_to_json(self, oracle_db_type):
        """ Converts oracle db type to an appropriate type for the JSON schema
        :param oracle_db_type: type in DGVa database
        :return: json_type: type in JSON schema
        """
        mapping_type = {
            "NUMBER(10,0)": "integer",
            "VARCHAR2(4000)": "string",
            #"NVARCHAR2": "string",
            "CLOB": "string",
            "DATE": "string",
            "TIMESTAMP": "string",
            "RAW": "string",
            "BLOB": "string",
            "BOOLEAN": "boolean",
            "NULL": "null",
        }
        json_type = mapping_type.get(oracle_db_type)
        return json_type

    def clean_data(self, df):
        """ Restricts data to the required values only. Changes type where needed.
        :param df: dataframe of metadata from DGVa
        :return: cleaned dataframe
        """
        clean_df = df.drop(columns=["Section"])
        clean_df_index = clean_df.set_index("metadataField")
        clean_df_index['type'] = clean_df_index['type'].map(self.map_type_oracle_to_json)
        fields_to_change = ["STUDY_UPDATE", "NEW_FEATURE"] # should be in upper case as not converted to camelCase
        matched_fields = clean_df_index.index.intersection(fields_to_change)
        clean_df_index.loc[matched_fields, "type"] = "boolean"
        return clean_df_index

    def generate_json_string(self, df):
        """Generates JSON object string. This string contains cleaned data originated from metadata spreadsheet.
        It will form the JSON object of a nested DGVA item.
        :param df: data frame
        :return: string of json object
        """
        json_string = df.to_json(orient="index", indent=4)
        return json_string

    def to_camel_case(self, string_to_convert):
        """ Convert string to camel case.
        :param string_to_convert: string
        :return camel case string
        """
        words = string_to_convert.replace("_", " ").replace("-", " ").split()
        if not words:
            return ""
        first_word = words[0].lower()
        other_words = [word.capitalize() for word in words[1:]]
        camel_case = first_word + "".join(other_words)
        return camel_case

    def nest_json_object(self, data_json_object):
        """ Turns a flat JSON object into a nested one.
        :param data_json_object: flat object
        :return: nested object, list of required fields
        """
        nested_object = {}
        required_fields = []
        for field_name, field_value in data_json_object.items():
            if field_value.get("notNull"):
                required_fields.append(field_name)
            camel_case_field_name = self.to_camel_case(field_name)
            nested_object[camel_case_field_name] = {
                "type": field_value.get("type"),
                "definition": field_value.get("definition"),
                "dbMapping": {
                    "dgvaTable": field_value.get("dgvaTable"),
                    "notNull": field_value.get("notNull"),
                    "columnDefault": field_value.get("default"),
                }
            }
        return nested_object, required_fields

    def create_json_schema(self):
        """This creates the JSON schema object.
        :return: JSON object for the schema
        """
        data_json_object = json.loads(self.data_string)
        nested_data_json_object, required_list = self.nest_json_object(data_json_object)


        schema = {
            "$schema": self.schema_version,
            "$id": self.id,
            "version": self.version,
            "author": self.author,
            "title": self.title,
            "description": self.description,
            "type": "object",
            "properties": {
                "dgva": {
                    "type": "array",
                    "description": "DGVa metadata that cannot be imported to EVA",
                    "items": {
                        "type": "object"
                        # --- this is where your object will be added ---
                    }
                }
            }
        }
        # adding the JSON object of a nested DGVA item to the schema
        schema["properties"]["dgva"]["items"]["properties"] = nested_data_json_object
        schema["properties"]["dgva"]["items"]["required"] = required_list
        return schema

    def print_json_schema(self, schema, output_path):
        """Prints the JSON schema to the output file.
        :param schema: JSON object
        :param output_path: string for path
        :return: None
        """
        with open(output_path, "w", newline="\n") as f:
            json.dump(schema, f, indent=4)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A tool to generate a JSON Schema for DGVa.")
    parser.add_argument("--input_path", help="Path to metadata tsv file")
    parser.add_argument("--output_path", help="Path to JSON schema file")
    args = parser.parse_args()

    path = args.input_path
    output_path = args.output_path
    builder = SchemaCreator(path)
    schema = builder.create_json_schema()
    builder.print_json_schema(schema, output_path)
