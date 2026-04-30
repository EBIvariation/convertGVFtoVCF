import os
import yaml


class ProjectPaths:
    """ The responsibility of this class is to manage paths from a config file.
    It resolves the full path for the files in the config.
    """
    def __init__(self, config_path="config.yaml"):
        cfg_path = config_path if config_path is not None else "config.yaml"
        self.root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) # TOP LEVEL
        self.package_dir = os.path.dirname(os.path.abspath(__file__)) # project directory of convertGVFtoVCF
        self.full_config_path = os.path.join(self.package_dir, "etc", cfg_path)
        # reading in
        try:
            with open(self.full_config_path, 'r') as f:
                data = yaml.safe_load(f) or {}
        except FileNotFoundError:
            print(f"Can't find file: {self.full_config_path}")
            data = {}
        paths = data.get('paths', {})
        self.etc_dir = os.path.join(self.package_dir, paths.get("etc_folder", ""))
        self.schema_path = os.path.join(self.package_dir, paths.get("schema_file", ""))
        self.test_dir = os.path.normpath(os.path.join(self.package_dir, paths.get("test_folder", "")))
