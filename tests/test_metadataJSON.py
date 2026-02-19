import json
from unittest import TestCase
from unittest.mock import patch, MagicMock, ANY

import oracledb

from convert_gvf_to_vcf.metadataJSON import DGVaMetadataRetriever

MOCK_CONFIG = {
    "DGVA": {
        "host": "mockhost",
        "port": "mockport",
        "user": "mockuser",
        "password": "mockpassword",
        "service_name": "mockservicename"
    }
}

class TestDGVaMetadataRetriever(TestCase):

    @patch('convert_gvf_to_vcf.metadataJSON.cfg.load_config_file')
    @patch('convert_gvf_to_vcf.metadataJSON.oracledb.connect', return_value=MagicMock())
    def test_connection(self, mock_connect, mock_config):
        mock_db = mock_config.return_value.DGVA.configure_mock(
            host=ANY, port=ANY, user=ANY, password=ANY, service_name=ANY
        )

        print(mock_db)
        client = DGVaMetadataRetriever(ANY)
        self.assertIsNotNone(client.connection)
        mock_connect.assert_called_once()

