import argparse

from convert_gvf_to_vcf.metadata_retrievers.dgvametadata import DGVAMetadataRetriever
from convert_gvf_to_vcf.metadata_retrievers.evametadata import EVAMetadataRetriever


def gather_metadata(config_input, json_output_eva, json_output_dgva, study_accession, vcf_output, assembly, assembly_report):

    retrieved_eva_metadata = EVAMetadataRetriever(config_input)
    with retrieved_eva_metadata:
        retrieved_eva_metadata.create_json_eva(json_file_path=json_output_eva, study_accession=study_accession, vcf_output=vcf_output, assembly=assembly, assembly_report=assembly_report)
    retrieved_dgva_metadata = DGVAMetadataRetriever(config_input)
    with retrieved_dgva_metadata:
        retrieved_dgva_metadata.create_json_dgva(json_output_dgva, study_accession)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config")
    parser.add_argument("--json_output_eva")
    parser.add_argument("--json_output_dgva")
    parser.add_argument("--study_accession")
    parser.add_argument("--vcf_output")
    parser.add_argument("--assembly")
    parser.add_argument("--assembly_report")
    args = parser.parse_args()

    gather_metadata(args.config, args.json_output_eva, args.json_output_dgva, args.study_accession, args.vcf_output, args.assembly, args.assembly_report)

if __name__ == "__main__":
    main()