import copy
import os
import json
from collections import defaultdict


from gather_metadata import gather_metadata_workflow, eva_update_metadata_with_vcf
from convertGVFtoVCF import convert, ProjectPaths
from ebi_eva_common_pyutils.logger import logging_config as log_cfg
logger = log_cfg.get_logger(__name__)

class GvfMetadataCoordinator:
    # the responsibility of this class is to co-ordinate multiple GVFs with their Metadata. Feeds into gather metadata script
    def __init__(self, study_and_gvf_files, base_output_dir, config_path):
        self.scan_results = study_and_gvf_files #{study_accession:[gvf_files]}
        self.base_output_dir = base_output_dir # top level output
        self.config_path = config_path
        self.project_paths = ProjectPaths()

    def process_studies(self):
        """ Processes the dictionary result of scanning the directories {study_accession:[gvf_files]}.
        If multiple files it will reconfigure the JSON file separated out by assembly.
        """
        empty_studies_log = []
        for study_accession, gvf_files in self.scan_results.items():
            # empty list of GVF files
            if not gvf_files:
                self._process_no_gvf_files(empty_studies_log, study_accession)
                continue
            # one file in list of GVFs
            if len(gvf_files) == 1:
                self._process_single_gvf_file(gvf_files, study_accession)
                continue
            # multiple gvf files
            elif len(gvf_files) > 1:
                self._process_multiple_gvf_files(gvf_files, study_accession)
                continue
            else:
                print("not in a recognised format")
        logger.info(f"Number of empty studies encountered: {len(empty_studies_log)}")

    def _process_multiple_gvf_files(self, gvf_files, study_accession):
        """Process if multiple gvf files are present.
        :params: gvf_files = list of gvf files
        :params: study_accession e.g estd1
        """
        logger.info(f"More than one GVF file found for {study_accession}. Separating by assembly.")
        # separate files by assembly {"GRCh37": [gvf_file_paths], "GRCh38": [gvf_file_paths]})
        assembly_groups = defaultdict(list)
        for gvf_file in gvf_files:
            file_name = os.path.basename(gvf_file)
            file_assembly = file_name.split(".")[2]
            assembly_groups[file_assembly].append(gvf_file)
        for assembly_name, files_in_assembly in assembly_groups.items():
            assembly_path, assembly_report_path, json_dgva, json_eva, vcf_output = self.set_up_inputs_and_outputs(
                assembly_name,
                study_accession)
            (eva_retriever, dgva_retriever, submitted_files, remapped_files) \
                = self.retrieve_metadata(json_eva, json_dgva, study_accession, assembly_path, assembly_report_path,
                                         files_in_assembly)
            # for those with multiple submitted files, reconfigure the JSON
            if len(submitted_files) > 1:
                self._determine_same_and_reconfigure_json(study_accession, submitted_files, json_eva, eva_retriever)

            # for those with multiple submitted files, reconfigure the JSON
            if len(remapped_files) > 1:
                self._determine_same_and_reconfigure_json(remapped_files, json_eva, eva_retriever)
            # TODO: check the resolution of GVF
            for individual_gvf in files_in_assembly:
                # Build unique VCF names based on the file name to avoid overwriting files
                base_name = os.path.basename(individual_gvf).replace(".gvf", "")
                individual_vcf_output = os.path.join(self.base_output_dir, f"{base_name}.vcf")
                convert(
                    gvf_input=individual_gvf,
                    vcf_output=individual_vcf_output,
                    assembly=assembly_path,
                    paths=self.project_paths
                )
                if eva_retriever:
                    eva_update_metadata_with_vcf(
                        eva_retriever=eva_retriever,
                        json_eva=json_eva,
                        vcf_output=individual_vcf_output
                    )

    def _process_single_gvf_file(self, gvf_files, study_accession):
        """Process if a single gvf files is present- get metadata as normal.
        :params: gvf_files = list of gvf files
        :params: study_accession e.g estd1
        """
        logger.info(f"One GVF file found for {study_accession}")
        gvf_file = gvf_files[0]
        assembly_name = gvf_file.split(".")[2]
        assembly_path, assembly_report_path, json_dgva, json_eva, vcf_output = self.set_up_inputs_and_outputs(
            assembly_name,
            study_accession)
        eva_retriever, dgva_retriever = self.retrieve_metadata(study_accession, assembly_name)
        # convert GVF to VCF
        convert(gvf_input=gvf_file, vcf_output=vcf_output, assembly=assembly_path, paths=self.project_paths)
        # add VCF details to the JSON file post-conversion
        if eva_retriever:
            eva_update_metadata_with_vcf(eva_retriever=eva_retriever, json_eva=json_eva, vcf_output=vcf_output)

    def _process_no_gvf_files(self, empty_studies_log, study_accession):
        """Process if no gvf files are present - logging.
        :params: empty_studies_log = list of empty studies
        :params: study_accession e.g estd1
        """
        logger.info(f"No GVF files found for {study_accession}")
        empty_studies_log.append(study_accession)

    def _determine_same_and_reconfigure_json(self, study_accession, target_files, json_eva, eva_retriever):
        """Determine if same biological or technical replicates, is so reconfigure JSON
        study_accession, target_files, json_eva, eva_retriever
        """
        determined_replicates = set()
        #TODO: explore gvf_file
        for gvf_file in target_files:
            # determine if they are the same (i.e. biological or technical replicates)
            reference = eva_retriever.referenceGenome or None
            analysis_types = eva_retriever._fetch_analysis_analysis_type(study_accession)
            method_types = eva_retriever._fetch_analysis_method_type(study_accession)
            experiment = eva_retriever._determine_analysis_experiment_type(analysis_types, method_types) or None

            replicate_signature = (reference, experiment, method_types)
            determined_replicates.add(replicate_signature)
        # If all files share the same metadata, reconfigure the JSON file in its analysis section
        if len(determined_replicates) == 1:
            self._reconfigure_json_multi_analysis(json_eva, target_files)
        else:
            # Leave json as it is (they differ so do not merge)
            pass

    def _reconfigure_json_multi_analysis(self, json_eva, gvf_files):
        """Reconfigures the EVA JSON with multiple files. This affects analysis, files and sample sections.
        :params: json_eva Path to EVA JSON file
        :params: files list of GVF files
        """
        try:
            with open(json_eva, 'r') as f_in:
                metadata = json.load(f_in)
            # find the important sections of the JSON file, leave the JSON file alone if not present.
            analysis_list = metadata.get("analysis", [])
            files_list = metadata.get("files", [])
            sample_list = metadata.get("sample", [])

            if not analysis_list or not files_list or not sample_list:
                return
            # update analysis and file blocks with new analysis aliases
            multiple_analyses, multiple_files, new_analysis_aliases = self._update_analysis_and_file_blocks(analysis_list, gvf_files,
                                                                                                            files_list)
            # after going through each GVF file, update the analysis alias in the sample block with the new analysis aliases
            self._update_sample_block(new_analysis_aliases, sample_list)

            # update the analysis, files and sample sections with the multiple datasets
            metadata["analysis"] = multiple_analyses
            metadata["files"] = multiple_files
            metadata["sample"] = sample_list

            # output to JSON
            with open(json_eva, 'w') as f_out:
                json.dump(metadata, f_out, indent=4)

        except (FileNotFoundError, json.JSONDecodeError) as err:
            print(f"Failed to update multi-analysis schema in JSON due to error: {err}")

    def _update_sample_block(self, new_analysis_aliases, sample_list):
        """Updates sample part of EVA JSON with new analysis aliases
        :params new_analysis_aliases: list of new names
        :params sample_list: block of EVA JSON
        """
        for sample_entry in sample_list:
            if "analysisAlias" in sample_entry and isinstance(sample_entry["analysisAlias"], list):
                sample_entry["analysisAlias"] = new_analysis_aliases

    def _update_analysis_and_file_blocks(self, analysis_list, files, files_list):
        """Updates analysis and file parts of the EVA JSON
        :params analysis_list: block of EVA JSON
        :params files: gvf giles
        :params files_list: block of EVA JSON
        :return multiple_analyses, multiple_files, new_analysis_alias: blocks for EVA JSON
        """
        # get the first blocks of the EVA JSON schema and expand on those
        initial_analysis_block = analysis_list[0]
        initial_file_block = files_list[0]

        multiple_analyses = []
        multiple_files = []
        new_analysis_aliases = []

        # change the Analysis and File blocks for each GVF file
        for index, file_path in enumerate(files, start=1):
            file_name = os.path.basename(file_path)

            # rename the analysis alias string
            base_alias = initial_analysis_block.get("analysisAlias", "analysis")
            unique_alias = f"{base_alias}_file_{index}"
            new_analysis_aliases.append(unique_alias)

            # copy the analysis block and update that with new analysis alias
            analysis_block = copy.deepcopy(initial_analysis_block)
            analysis_block["analysisAlias"] = unique_alias
            multiple_analyses.append(analysis_block)

            # copy the file block and update that with new analysis alias
            file_block = copy.deepcopy(initial_file_block)
            file_block["analysisAlias"] = unique_alias
            file_block["fileName"] = file_name
            multiple_files.append(file_block)
        return multiple_analyses, multiple_files, new_analysis_aliases

    def retrieve_metadata(self, json_eva, json_dgva, study_accession, assembly_path, assembly_report_path, files_in_assembly = None):
        """Retrieves metadata and if applicable, returns list of submitted and remapped files
        :params  json_eva, json_dgva, study_accession, assembly_path, assembly_report_path, files_in_assembly = None
        :returns eva_retriever, dgva_retriever or eva_retriever, dgva_retriever, submitted_files, remapped_files
        """
        # fetch the metadata
        eva_retriever, dgva_retriever = gather_metadata_workflow(
            config=self.config_path,
            json_eva=json_eva,
            json_dgva=json_dgva,
            study_accession=study_accession,
            assembly=assembly_path,
            assembly_report=assembly_report_path
        )
        if files_in_assembly:
            # for each assembly version, check if multiple submitted files
            submitted_files = [f for f in files_in_assembly if "Submitted" in os.path.basename(f)]
            remapped_files = [f for f in files_in_assembly if "Remapped" in os.path.basename(f)]
            return eva_retriever, dgva_retriever, submitted_files, remapped_files
        else:
            return eva_retriever, dgva_retriever

    def set_up_inputs_and_outputs(self, assembly_name, study_accession):
        """Programmatically sets up input and output paths.
        :params assembly_name e.g. GRCh37
        :params study_accession e.g. estd1
        :return assembly_path, assembly_report_path, json_dgva, json_eva, vcf_output: input, input, output, output, output
        """
        assembly_name_formatted = f"_{assembly_name}_" if assembly_name else "_"
        # get the inputs
        assembly_path = self.project_paths.assembly_paths.get(assembly_name)
        assembly_report_path = self.project_paths.assembly_report_paths.get(assembly_name)

        # programmatically generate the output files by assembly
        json_eva = os.path.join(self.base_output_dir, f"{study_accession}{assembly_name_formatted}eva.json")
        json_dgva = os.path.join(self.base_output_dir, f"{study_accession}{assembly_name_formatted}dgva.json")
        vcf_output = os.path.join(self.base_output_dir, f"{study_accession}{assembly_name_formatted}.vcf")
        return assembly_path, assembly_report_path, json_dgva, json_eva, vcf_output

