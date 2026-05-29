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



        def scan(self):
            """Scans entire HPC.
            :returns: dictionary of lists {study_accession:[list_of_gvf_files]}
            """
            logger.info(f"Scanning HPC directory: {self.hpc_dir}")
            all_extensions = set()
            study_and_files = dict()
            for (current_dir, dirs, files) in os.walk(self.hpc_dir):
                # remove hidden files
                dirs[:] = [d for d in dirs if not d.startswith('.')]
                folder_name = os.path.basename(current_dir)
                if current_dir == self.hpc_dir:
                    logger.info(f"Study accessions found in HPC directory: {len(dirs)}")
                elif folder_name.startswith(("estd", "nstd")):
                    study_accession_of_folder = folder_name.split("_")[0]
                    if len(dirs) > 1:
                        logger.warning(f"Study accession {study_accession_of_folder} has more than 1 directory")
                elif folder_name == "gvf":
                    logger.info(f"Found {study_accession_of_folder}: {len(files)} files.")
                    for file in files:
                        logger.info(f"File: {file}")
                        ext = os.path.splitext(file)[1].lower()
                        all_extensions.add(ext)
                    if files and all_extensions != {".gvf"}:
                        logger.warning(f"Files in study accession {study_accession_of_folder} are not all GVFs.")
                        gvf_files = [file for file in files if file.endswith(".gvf")]
                        study_and_files[study_accession_of_folder] = gvf_files
                    else:
                        study_and_files[study_accession_of_folder] = files
                else:
                    logger.warning(f"Unconventional directory structure: {current_dir}.")
            return study_and_files

        def find_files(self):
            pass

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
    gvf_data = finder.scan()
    print(gvf_data)
    print("complete")


if __name__ == "__main__":
    main()