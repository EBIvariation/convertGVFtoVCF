import argparse
import json
import os.path
import shutil

from convert_gvf_to_vcf.metadata_retrievers.dgvametadata import DGVAMetadataRetriever
from convert_gvf_to_vcf.metadata_retrievers.evametadata import EVAMetadataRetriever
from ebi_eva_common_pyutils.logger import logging_config as log_cfg
logger = log_cfg.get_logger(__name__)

def eva_gather_metadata(retriever, json_output, study_accession, assembly, assembly_report):
    retriever.create_json_eva(json_file_path=json_output, study_accession=study_accession, assembly=assembly, assembly_report=assembly_report)

def eva_add_file_metadata(retriever, json_output, vcf_output):
    files_file_name = retriever._get_file_name(vcf_output)
    files_file_size = retriever._get_file_size(vcf_output)
    files_file_md5 = retriever._get_file_md5(vcf_output)

    try:
        with open(json_output, 'r') as f_in:
            metadata = json.load(f_in)
    except FileNotFoundError:
        logger.error(f"Cannot find json output: {json_output}")
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
    parser.add_argument("--json_output_eva")
    parser.add_argument("--json_output_dgva")
    parser.add_argument("--study_accession")
    parser.add_argument("--vcf_output")
    parser.add_argument("--assembly")
    parser.add_argument("--assembly_report")
    args = parser.parse_args()

    eva_retriever = EVAMetadataRetriever(
        path_to_config_yaml=args.config
    )
    eva_gather_metadata(eva_retriever, args.json_output_eva, args.study_accession, args.assembly, args.assembly_report)
    eva_add_file_metadata(eva_retriever, args.json_output_eva, args.vcf_output)

    dgva_retriever = DGVAMetadataRetriever(
        path_to_config_yaml=args.config
    )
    dgva_retriever.create_json_dgva(json_file_path=args.json_output_dgva, study_accession=args.study_accession)



if __name__ == "__main__":
    main()