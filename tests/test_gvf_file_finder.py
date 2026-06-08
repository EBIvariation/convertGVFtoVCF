import unittest
import os
from unittest.mock import MagicMock, patch

from convert_gvf_to_vcf.gvf_file_finder import GvfFileFinder
from convert_gvf_to_vcf.project_paths import ProjectPaths


class TestGvfFileFinder(unittest.TestCase):
    def setUp(self):
        self.paths = ProjectPaths()
        self.top_dir = os.path.join(self.paths.test_dir, "data_dir")
        self.current_dir = os.path.join(self.top_dir, "estd1_Redon_et_al_2006", "gvf")
        self.study_accession = "estd1"

    def test_process_gvf_directory(self):
        all_extensions = set()
        current_dir = self.current_dir
        dirs = list()
        files = ['estd1_Redon_et_al_2006.2014-04-01.GRCh37.Remapped.gvf', 'estd1_Redon_et_al_2006.2014-04-01.NCBI35.Submitted.gvf', 'estd1_Redon_et_al_2006.2014-04-01.GRCh37.p13.Remapped.gvf', 'estd1_Redon_et_al_2006.2014-04-01.GRCh38.Remapped.gvf']
        study_and_files = dict()
        study_file_dict = GvfFileFinder._process_gvf_directory(self, all_extensions=all_extensions, current_dir=current_dir, dirs=dirs, files=files, study_and_files=study_and_files, study_accession=None)
        assert len(study_file_dict) == 1
        assert self.study_accession in study_file_dict
        expected_paths = []
        for file in files:
            expected_paths.append(os.path.join(self.current_dir, file))
        assert study_file_dict[self.study_accession] == expected_paths

    def test_deduplicate_files(self):
        files = ['estd1_Redon_et_al_2006.2014-04-01.GRCh37.Remapped.gvf', 'estd1_Redon_et_al_2006.2014-04-01.NCBI35.Submitted.gvf', 'estd1_Redon_et_al_2006.2014-04-01.GRCh37.p13.Remapped.gvf', 'estd1_Redon_et_al_2006.2014-04-01.GRCh38.Remapped.gvf']
        expected_paths = []
        for file in files:
            expected_paths.append(os.path.join(self.current_dir, file))
        study_and_files = {self.study_accession: expected_paths}
        self._check_md5_match = MagicMock(return_value=False)
        os.path.getsize = MagicMock(return_value=100)
        result = GvfFileFinder.deduplicate_files(self, study_and_files)
        assert "estd1" in result
        assert result["estd1"] == expected_paths

    @patch('os.walk')
    def test_scan(self, mock_walk):
        file_finder = GvfFileFinder(search_dir=self.top_dir)
        study_name = "estd1_Redon_et_al_2006"

        study_path = os.path.join(self.top_dir, study_name)
        gvf_path = os.path.join(study_path, "gvf")

        files = ['estd1_Redon_et_al_2006.2014-04-01.GRCh37.Remapped.gvf', 'estd1_Redon_et_al_2006.2014-04-01.NCBI35.Submitted.gvf', 'estd1_Redon_et_al_2006.2014-04-01.GRCh37.p13.Remapped.gvf', 'estd1_Redon_et_al_2006.2014-04-01.GRCh38.Remapped.gvf']
        expected_paths = []
        for file in files:
            expected_paths.append(os.path.join(gvf_path, file))
        mock_return_dict = {self.study_accession: expected_paths} # assuming no duplicates paths
        file_finder.deduplicate_files = MagicMock(return_value=mock_return_dict)
        file_finder._process_gvf_directory = MagicMock(return_value =mock_return_dict)

        mock_walk.return_value = [
            (self.top_dir, [study_name], []),
            (study_path, ["gvf"], []),
            (gvf_path, [], files)
        ]
        result = file_finder.scan(self.study_accession)
        self.assertEqual(result[self.study_accession], expected_paths)
        file_finder._process_gvf_directory.assert_called_once()
        file_finder.deduplicate_files.assert_called_once()

    def test_get_md5(self):
        file_to_test = os.path.join(self.current_dir, 'estd1_Redon_et_al_2006.2014-04-01.GRCh37.Remapped.gvf')
        actual_md5 = GvfFileFinder.get_md5(self, file_to_test)
        expected_md5 = "eccbc87e4b5ce2fe28308fd9f2a7baf3"
        self.assertEqual(actual_md5, expected_md5)