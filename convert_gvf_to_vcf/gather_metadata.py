import argparse
import json
import os.path
import shutil

from convert_gvf_to_vcf.metadataJSON import DGVaMetadataRetriever


def gather_metadata(retriever, json_output, study_accession, assembly, assembly_report):
    retriever.create_json_file(json_file_path=json_output, study_accession=study_accession, assembly=assembly, assembly_report=assembly_report)

def add_file_metadata(retriever, json_output,vcf_output):
    files_file_name = retriever._get_file_name(vcf_output)
    files_file_size = retriever._get_file_size(vcf_output)
    files_file_md5 = retriever._get_file_md5(vcf_output)

    with open(json_output, 'r') as f_in:
        metadata = json.load(f_in)

    # moving the file to a preconversion file to prevent confusion
    base_path, ext = os.path.splitext(json_output)
    preconversion_json_path = base_path + "_preconverted.json" # this will be missing part of the files section
    shutil.copy(json_output, preconversion_json_path)

    # adding the missing files sections
    for file_object in metadata["files"]:
        file_object["fileName"] = files_file_name
        file_object["fileSize"] = files_file_size
        file_object["md5"] = files_file_md5

    with open(json_output, 'w') as f_out:
        json.dump(metadata, f_out, indent=4)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config")
    parser.add_argument("--path_config")
    parser.add_argument("--json_output")
    parser.add_argument("--study_accession")
    parser.add_argument("--vcf_output")
    parser.add_argument("--assembly")
    parser.add_argument("--assembly_report")
    args = parser.parse_args()

    retriever = DGVaMetadataRetriever(
        path_to_config_yaml=args.config,
        path_to_path_config_yaml=args.path_config
    )
    gather_metadata(retriever, args.json_output, args.study_accession, args.assembly, args.assembly_report)
    add_file_metadata(retriever, args.json_output, args.vcf_output)
if __name__ == "__main__":
    main()