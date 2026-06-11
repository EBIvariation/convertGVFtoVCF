import os.path
import unittest

from convert_gvf_to_vcf.project_paths import ProjectPaths


class TestGvfMetadataCoordinator(unittest.TestCase):
    def setUp(self):
        self.paths = ProjectPaths()
        self.config = self.paths.full_config_path
        self.test_dir = self.paths.test_dir
        self.output = os.path.join(self.test_dir, "output", "gvf_metadata_coord")



if __name__ == '__main__':
    unittest.main()
