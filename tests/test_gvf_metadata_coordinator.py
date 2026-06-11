import os.path
import unittest
from unittest.mock import patch, MagicMock

from convert_gvf_to_vcf.project_paths import ProjectPaths
from convert_gvf_to_vcf.gather_metadata import (
    eva_update_metadata_with_vcf,
    gather_metadata_workflow,
)
from convert_gvf_to_vcf.gvf_metadata_coordinator import GvfMetadataCoordinator


class TestGvfMetadataCoordinator(unittest.TestCase):
    def setUp(self):
        self.paths = ProjectPaths()
        self.config = self.paths.full_config_path
        self.test_dir = self.paths.test_dir
        self.output_dir = os.path.join(self.test_dir, "output", "gvf_metadata_coord")


    @patch('convert_gvf_to_vcf.gvf_metadata_coordinator.logger')
    def test_process_studies(self, mock_logger):
        # testing no GVF files
        study_and_gvf_files_input_data = {"estd1": []}
        coordinator = GvfMetadataCoordinator(study_and_gvf_files_input_data, self.output_dir, self.config)
        coordinator._process_no_gvf_files = MagicMock()
        coordinator.process_studies()
        coordinator._process_no_gvf_files.assert_called_once_with([], "estd1")
        # testing one GVF file
        study_and_gvf_files_input_data = {"estd1": ["path/to/file1.gvf"]}
        coordinator = GvfMetadataCoordinator(study_and_gvf_files_input_data, self.output_dir, self.config)
        coordinator._process_single_gvf_file = MagicMock()
        coordinator.process_studies()
        coordinator._process_single_gvf_file.assert_called_once_with(["path/to/file1.gvf"], "estd1")
        # testing multiple GVF file
        study_and_gvf_files_input_data = {"estd1": ["path/to/file1.gvf", "path/to/file2.gvf"]}
        coordinator = GvfMetadataCoordinator(study_and_gvf_files_input_data, self.output_dir, self.config)
        coordinator._process_multiple_gvf_files = MagicMock()
        coordinator.process_studies()
        coordinator._process_multiple_gvf_files.assert_called_once_with(["path/to/file1.gvf", "path/to/file2.gvf"], "estd1")

if __name__ == '__main__':
    unittest.main()
