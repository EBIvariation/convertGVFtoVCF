import argparse

from convert_gvf_to_vcf.metadataJSON import DGVaMetadataRetriever


def gather_metadata(config_input, json_output, study_accession, vcf_output, assembly, assembly_report):

    retrieved_dgva_metadata = DGVaMetadataRetriever(config_input)
    with retrieved_dgva_metadata:
        retrieved_dgva_metadata.create_json_file(json_file_path=json_output, study_accession=study_accession, vcf_output=vcf_output, assembly=assembly, assembly_report=assembly_report)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config")
    parser.add_argument("--json_output")
    parser.add_argument("--study_accession")
    parser.add_argument("--vcf_output")
    parser.add_argument("--assembly")
    parser.add_argument("--assembly_report")
    args = parser.parse_args()

    gather_metadata(args.config, args.json_output, args.study_accession, args.vcf_output, args.assembly, args.assembly_report)

if __name__ == "__main__":
    main()