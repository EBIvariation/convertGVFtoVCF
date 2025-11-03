"""This sets up the logger for logging messages"""
import logging
import os

def set_up_logging(log_path=None):
    # Create logger for module name
    logger = logging.getLogger(__name__)

    # Set up paths
    input_folder = os.path.dirname(__file__)
    output_folder = os.path.join(input_folder, "..","tests", "output")
    if log_path is None:
        log_path = os.path.join(output_folder, "converted.log")
    logging.basicConfig(filename=log_path,
                        level=logging.DEBUG,
                        format="%(asctime)s - %(levelname)s - %(message)s")
    return logger, log_path