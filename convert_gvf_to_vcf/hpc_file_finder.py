import argparse
import os
from ebi_eva_common_pyutils.config import cfg
from ebi_eva_common_pyutils.logger import logging_config as log_cfg

from convert_gvf_to_vcf.project_paths import ProjectPaths
from convert_gvf_to_vcf.utils import get_validated_value

logger = log_cfg.get_logger(__name__)
class HpcFileFinder:
        """ The responsibility of this class is to traverse the FTP directory to find the correct GVF file paths."""
        def __init__(self, ftp_dir, study_accession):
            self.paths = ProjectPaths()
            self.hpc_dir = ftp_dir
            self.study_accession = study_accession

        def _extract_target_tuple(self, study_accession):
            """Creates a target tuple using a study accession
            :params: study_accession
            :return: target_study_accession (tuple)
            """
            if study_accession:
                if isinstance(study_accession, str):
                    target_study_accession = (study_accession,)
                else:
                    target_study_accession = tuple(study_accession)
            else:
                target_study_accession = None
            return target_study_accession

        def _process_gvf_directory(self, all_extensions, current_dir, dirs, files, study_and_files,
                                   target_study_accession):
            """ Processes the GVF directory (e.g. /hpc_dir/estd3_Name_et_al_2008/gvf). Returns the study accession and the GVF files.
            :params all_extensions: set of file extensions
            :params current_dir: from os.walk - root
            :params dirs: from os.walk - subdirectories
            :params files: from os.walk - files
            :params study_and_files: dictionary
            :params target_study_accessions: the study accession restricted to if applicable
            :return study_and_files: dictionary of lists {study_accession:[list_of_gvf_full_paths]}
            """
            parent_dir = os.path.basename(os.path.dirname(current_dir))
            study_accession_of_folder = parent_dir.split("_")[0]
            if target_study_accession:
                logger.info(
                    f"AFTER restricting target = {len(dirs)} study accessions found in HPC directory: {study_accession_of_folder}")
            logger.info(f"Found {study_accession_of_folder}: {len(files)} files.")
            for file in files:
                logger.info(f"File: {file}")
                ext = os.path.splitext(file)[1].lower()
                all_extensions.add(ext)
            if files and all_extensions != {".gvf"}:
                logger.warning(f"Files in study accession {study_accession_of_folder} are not all GVFs.")
                gvf_files = [os.path.join(current_dir, file) for file in files if file.endswith(".gvf")]
                study_and_files[study_accession_of_folder] = gvf_files
            else:
                study_and_files[study_accession_of_folder] = [os.path.join(current_dir, file) for file in files]
            return study_and_files

        def scan(self, study_accession=None):
            """Scans entire HPC.
            :returns: dictionary of lists {study_accession:[list_of_gvf_files]}
            """
            logger.info(f"Scanning HPC directory: {self.hpc_dir}")
            study_and_files = dict()
            # set target tuple
            target_study_accession = self._extract_target_tuple(study_accession)
            for (current_dir, dirs, files) in os.walk(self.hpc_dir):
                # remove hidden files
                dirs[:] = [d for d in dirs if not d.startswith('.')]
                # go through hpc dir of /hpc_dir/estd3_Name_et_al_2008/gvf
                if current_dir == self.hpc_dir:
                    if target_study_accession is None:
                        logger.info(f"Study accessions found in HPC directory: {len(dirs)}")
                    # set restriction
                    if target_study_accession:
                        logger.info(f"BEFORE restricting target = {len(dirs)} Study accessions found in HPC directory.")
                        dirs[:] = [d for d in dirs if d.startswith(target_study_accession)]
                        logger.info(f"Target study accession set to {target_study_accession}")
                    continue
                all_extensions = set()
                folder_name = os.path.basename(current_dir)
                # go through estd3_Name_et_al_2008 of /hpc_dir/estd3_Name_et_al_2008/gvf
                if folder_name.startswith(("estd", "nstd")):
                    study_accession_of_folder = folder_name.split("_")[0]
                    if len(dirs) > 1:
                        logger.warning(f"Study accession {study_accession_of_folder} has more than 1 directory")
                # go through gvf of /hpc_dir/estd3_Name_et_al_2008/gvf
                elif folder_name == "gvf":
                    study_and_files = self._process_gvf_directory(all_extensions, current_dir, dirs, files, study_and_files,
                                                target_study_accession)
                else:
                    logger.warning(f"Unconventional directory structure: {current_dir}.")
            return study_and_files




def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--ftp_dir", required=True, help="FTP dir.")
    parser.add_argument("--log", required=True, help="Path to log")
    parser.add_argument("--study_accession", help="Study accession.")

    args = parser.parse_args()

    # Set up logging functionality
    if args.log:
        log_cfg.add_file_handler(args.log)
        logger.info(f"The log file is {args.log}")
    else:
        log_cfg.add_stdout_handler()

    finder = HpcFileFinder(ftp_dir = args.ftp_dir, study_accession= args.study_accession)
    gvf_data = finder.scan(args.study_accession)
    print(gvf_data)
    print("complete")


if __name__ == "__main__":
    main()