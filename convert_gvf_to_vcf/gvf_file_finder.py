import argparse
import hashlib
import os


from ebi_eva_common_pyutils.config import cfg
from ebi_eva_common_pyutils.logger import logging_config as log_cfg

from convert_gvf_to_vcf.project_paths import ProjectPaths
from convert_gvf_to_vcf.utils import get_validated_value

logger = log_cfg.get_logger(__name__)
class GvfFileFinder:
        """ The responsibility of this class is to traverse the directory to find the correct GVF file paths."""
        def __init__(self, search_dir):
            self.paths = ProjectPaths()
            self.search_dir = search_dir

        def _process_gvf_directory(self, all_extensions, current_dir, dirs, files, study_and_files,
                                   study_accession):
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
            if study_accession:
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
            :returns: dictionary of lists {study_accession:[unique_gvf_files]}
            """
            logger.info(f"Scanning directory: {self.search_dir}")
            study_and_files = dict()
            # set target tuple
            for (current_dir, dirs, files) in os.walk(self.search_dir):
                # remove hidden files
                dirs[:] = [d for d in dirs if not d.startswith('.')]
                # go through hpc dir of /hpc_dir/estd3_Name_et_al_2008/gvf
                if current_dir == self.search_dir:
                    if study_accession is None:
                        logger.info(f"Study accessions found in directory: {len(dirs)}")
                    # set restriction
                    if study_accession:
                        logger.info(f"BEFORE restricting target = {len(dirs)} Study accessions found in HPC directory.")
                        dirs[:] = [d for d in dirs if d.startswith(study_accession)]
                        logger.info(f"Target study accession set to {study_accession}")
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
                                                study_accession)
                else:
                    logger.warning(f"Unconventional directory structure: {current_dir}.")
            study_and_gvf_files = self.deduplicate_files(study_and_files)
            return study_and_gvf_files

        def deduplicate_files(self, study_and_files):
            """ Deduplicate the GVF files found in the dictionary of lists {study_accession: [listofGVFs]}
            :params study_and_files: dictionary of list with duplicated elements
            :returns study_and_unique_files: dictionary of list with unique elements
            """
            study_and_unique_files = {}
            for study_accession, gvf_list in study_and_files.items():
                seen_files = {}
                unique_files = []
                for gvf in gvf_list:
                    file_size = os.path.getsize(gvf)
                    files_with_the_same_size = seen_files.get(file_size, [])
                    matched_file = next((existing_file for existing_file in files_with_the_same_size if self._check_md5_match(gvf, existing_file)), None)
                    if matched_file:
                        logger.error(f"Duplicate files have been found: {gvf} and {matched_file}")
                    else:
                        unique_files.append(gvf)
                        seen_files[file_size] = files_with_the_same_size + [gvf]
                study_and_unique_files[study_accession] = unique_files
            return study_and_unique_files

        def get_md5(self, file_path):
            """Gets md5 of a file path
            :params: file_path: gvf file
            :returns: md5
            """
            hasher = hashlib.md5()
            try:
                with open(file_path, "rb") as f:
                    for chunk in iter(lambda: f.read(65536), b""):
                        hasher.update(chunk)
                return hasher.hexdigest()
            except (OSError, IOError):
                return None

        def _check_md5_match(self, file_path_1, file_path_2):
            """Returns True if both files have the exact same MD5 hash."""
            md5_path1 = self.get_md5(file_path_1)
            md5_path2 = self.get_md5(file_path_2)
            if md5_path1 and md5_path2:
                return self.get_md5(file_path_1) == self.get_md5(file_path_2)
            return False



def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--search_dir", required=True, help="Directory to search i.e. FTP dir.")
    parser.add_argument("--log", required=True, help="Path to log")
    parser.add_argument("--study_accession", help="Study accession.")

    args = parser.parse_args()

    # Set up logging functionality
    if args.log:
        log_cfg.add_file_handler(args.log)
        logger.info(f"The log file is {args.log}")

    finder = GvfFileFinder(search_dir= args.search_dir)
    gvf_data = finder.scan(args.study_accession)
    print(gvf_data)
    print("complete")


if __name__ == "__main__":
    main()