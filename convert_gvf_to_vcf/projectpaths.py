import os
import yaml


class ProjectPaths:
    """ The responsibility of this class is to manage paths from a config file.
    It resolves the full path for the files in the config.
    """
    def __init__(self, config_path="config.yaml"):
        self.base_dir = os.path.dirname(os.path.abspath(__file__)) # project directory of convertGVFtoVCF
        self.full_config_path = os.path.join(self.base_dir, "etc", config_path)
        with open(self.full_config_path, 'r') as f:
            data = yaml.safe_load(f)
        self.etc_dir = os.path.join(self.base_dir, data['paths']['etc_folder'])
        self.schema_path = os.path.join(self.base_dir, data['paths']['schema_file'])

