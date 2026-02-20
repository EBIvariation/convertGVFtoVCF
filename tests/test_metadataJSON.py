import json
import os
from unittest import TestCase
from unittest.mock import patch, MagicMock, ANY, Mock

import oracledb
from ebi_eva_common_pyutils.spreadsheet.metadata_xlsx_utils import metadata_xlsx_version

from convert_gvf_to_vcf.metadataJSON import DGVaMetadataRetriever


class TestDGVaMetadataRetriever(TestCase):
    def setUp(self):
        input_folder = os.path.dirname(__file__)
        self.config = os.path.join(input_folder, "input", "test.config")

    def test_init(self):
        # testing the __init__ function has loaded the config correctly
        metadata_client = DGVaMetadataRetriever(self.config)
        assert metadata_client._host == "mockhost"
        assert metadata_client._port == "1111"
        assert metadata_client._username == "mockuser"
        assert metadata_client._password == "mockpassword"
        assert metadata_client._service_name == "mockservicename"
        assert metadata_client._connection == None
        assert metadata_client._max_retries == 3

    # TARGET OBJECT = convert_gvf_to_vcf.metadataJSON.oracledb.connect
    # PATCH = @patch
    # MOCK OBJECT = mock_connection
    @patch("convert_gvf_to_vcf.metadataJSON.oracledb.connect")
    def test_connection_new(self, mock_connection):
        mock_connection_object = Mock()
        mock_connection.return_value = mock_connection_object
        metadata_client = DGVaMetadataRetriever(self.config)
        result = metadata_client.connection
        self.assertEqual(result, mock_connection_object)
        mock_connection.assert_called_once()

    # TARGET OBJECT = convert_gvf_to_vcf.metadataJSON.oracledb.connect
    # PATCH = @patch
    # MOCK OBJECT = mock_connection
    @patch("convert_gvf_to_vcf.metadataJSON.oracledb.connect")
    def test_connection_healthy(self, mock_connection):
        # expected behaviour: if healthy, keep the connection
        mock_connection_object_existing = Mock()
        mock_connection.return_value = mock_connection_object_existing
        metadata_client = DGVaMetadataRetriever(self.config)
        metadata_client._connection = mock_connection_object_existing
        result = metadata_client.connection
        self.assertEqual(result, mock_connection_object_existing)
        mock_connection.assert_not_called()

    # TARGET OBJECT = convert_gvf_to_vcf.metadataJSON.oracledb.connect
    # PATCH = @patch
    # MOCK OBJECT = mock_connection
    @patch("convert_gvf_to_vcf.metadataJSON.oracledb.connect")
    def test_connection_unhealthy(self, mock_connection):
        # expected behaviour: if unhealthy, close the connection then open a new connection
        mock_unhealthy = Mock()
        mock_unhealthy.is_healthy.return_value = False
        mock_new_connection = Mock()
        mock_connection.return_value = mock_new_connection

        metadata_client = DGVaMetadataRetriever(self.config)
        # set connection to unhealthy
        metadata_client._connection = mock_unhealthy
        result = metadata_client.connection
        # did you close the unhealthy connection
        mock_unhealthy.close_assert_called_once()
        # did you get a new connection
        self.assertEqual(result, mock_new_connection)

    # TARGET OBJECT = convert_gvf_to_vcf.metadataJSON.oracledb.connect
    # PATCH = @patch
    # MOCK OBJECT = mock_connection
    @patch("convert_gvf_to_vcf.metadataJSON.oracledb.connect")
    def test_connection_max_retry_third_attempt_success(self, mock_connection):
        # expected behaviour: unsuccessful on first two attempts, successful on third attempt.
        mock_healthy = Mock()
        mock_healthy.is_healthy.return_value = True
        mock_connection.side_effect = [Exception, Exception, mock_healthy]

        metadata_client = DGVaMetadataRetriever(self.config)
        # checking attribute is set to max 3 retries
        self.assertEqual(metadata_client._max_retries, 3)
        result = metadata_client.connection
        # checking your connection is healthy
        self.assertEqual(result, mock_healthy)
        self.assertEqual(mock_connection.call_count, 3)

    # TARGET OBJECT = convert_gvf_to_vcf.metadataJSON.oracledb.connect
    # PATCH = @patch
    # MOCK OBJECT = mock_connection
    @patch("convert_gvf_to_vcf.metadataJSON.oracledb.connect")
    def test_connection_max_retry_third_attempt_fail(self, mock_connection):
        # expected behaviour: it will allow 3 retries. it will not allow the fourth attempt
        mock_connection.side_effect = [Exception("attempt1"), Exception("attempt2"), Exception("attempt3"), Exception("attempt4")]
        metadata_client = DGVaMetadataRetriever(self.config)
        with self.assertRaises(Exception) as unsuccessful:
            result = metadata_client.connection
        # this should only allow 3 attempts
        assert mock_connection.call_count == 3
