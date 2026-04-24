import os
import yaml


class ProjectPaths:
    """ The responsibility of this class is to manage paths from a config file.
    It resolves the full path for the files in the config.
    """
    def __init__(self, config_path="config.yaml"):
        self.root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) # TOP LEVEL
        self.package_dir = os.path.dirname(os.path.abspath(__file__)) # project directory of convertGVFtoVCF

        self.full_config_path = os.path.join(self.package_dir, "etc", config_path)
        # reading in
        with open(self.full_config_path, 'r') as f:
            data = yaml.safe_load(f)
        self.etc_dir = os.path.join(self.package_dir, data['paths']['etc_folder'])
        self.schema_path = os.path.join(self.package_dir, data['paths']['schema_file'])
        self.test_dir = os.path.normpath(os.path.join(self.package_dir, data['paths']['test_folder']))
        print("base_dir", self.package_dir)
ProjectPaths()