import json
import os
from unittest import TestCase
from unittest.mock import patch, MagicMock, ANY, Mock
from collections import namedtuple
import oracledb
from ebi_eva_common_pyutils.spreadsheet.metadata_xlsx_utils import metadata_xlsx_version

from convert_gvf_to_vcf.metadataJSON import DGVaMetadataRetriever


class TestDGVaMetadataRetriever(TestCase):
    def setUp(self):
        input_folder = os.path.dirname(__file__)
        self.config = os.path.join(input_folder, "input", "test.config")
        self.vcf_output = os.path.join(input_folder, "output", "a.vcf")

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
        mock_unhealthy.close.assert_called_once()
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

    def test_load_from_db(self):
        pass

    def test_create_json_file(self):
        pass

    def test__get_validated_value(self):
        pass


    @patch("convert_gvf_to_vcf.metadataJSON.DGVaMetadataRetriever._fetch_submitter_details_dictionary")
    def test_get_submitter_details(self, mock_fetch):

        # first_name, last_name, telephone, email, laboratory, centre, address
        # create mock data from load_from_db
        mock_data = {
            "0": ["John", "Smith", "0123456789", "e1@mail.com",  "centre1", "123 Road Name, City, Country"],  # Test case: positive test case, expected inputs
            "1": ["Jane", None, "0123456789", "e2@mail.com",  "centre2","123 Road Name, City, Country"],      # Test Case: Last Name - None
            "2": ["Bob", "", "0123456789", "e3@mail.com" , "centre3", "123 Road Name, City, Country"],       # Test Case: Last Name - empty string
            "3": ["Mike", "Jones", "0123456789", None ,  "centre4", "123 Road Name, City, Country"],          # Test Case: Email Add - None
            "4": ["Dave", "Smith", "0123456789", "" , "centre4", "123 Road Name, City, Country"],            # Test Case: Email Add - empty string
            "5": ["Bobby", "Jones", "0123456789", "e4@mail.com" ,  None, "123 Road Name, City, Country"],     # Test Case: Centre - None
            "6": ["Mary", "Smith", None, "e1@mail.com", "centre1", "123 Road Name, City, Country"],          # Test Case Telephone - None
            "7": ["Robert", "Williams", "", "e1@mail.com", "centre1", "123 Road Name, City, Country"],       # Test Case Telephone - empty string
            "8": ["Patricia", "Smith", "", "e1@mail.com", "centre1", None],                                  # Test Case Address - None
            "9": ["John", "Brown", "", "e1@mail.com", "centre1", ""]                                         # Test Case Address - empty string
        }
        mock_fetch.return_value = mock_data

        metadata_client = DGVaMetadataRetriever(self.config)
        result = metadata_client._get_submitter_details("STUDY123")

        # 3. ASSERT: Verify the default values were applied
        expected = [{
                'firstName': 'John',
                'lastName': 'Smith',
                'telephone': '0123456789',
                'laboratory': '',
                'email': 'e1@mail.com',
                'centre': 'centre1',
                'address': '123 Road Name, City, Country'
            },
            {
                'firstName': 'Jane',
                'lastName': '',
                'telephone': '0123456789',
                'laboratory': '',
                'email': 'e2@mail.com',
                'centre': 'centre2',
                'address': '123 Road Name, City, Country'
             },
            {
                'firstName': 'Bob',
                'lastName': '',
                'telephone': '0123456789',
                'laboratory': '',
                'email': 'e3@mail.com',
                'centre': 'centre3',
                'address': '123 Road Name, City, Country'
            },
            {
                'firstName': 'Mike',
                'lastName': 'Jones',
                'telephone': '0123456789',
                'centre': 'centre4',
                'email': '',
                'laboratory': '',
                'address': '123 Road Name, City, Country'
            },
            {
                'firstName': 'Dave',
                'lastName': 'Smith',
                'telephone': '0123456789',
                'email': '',
                'laboratory': '',
                'centre': 'centre4',
                'address': '123 Road Name, City, Country'
            },
            {
                'firstName': 'Bobby',
                'lastName': 'Jones',
                'telephone': '0123456789',
                'email': 'e4@mail.com',
                'laboratory': '',
                'centre': '',
                'address': '123 Road Name, City, Country'
            },
            {
                'address': '123 Road Name, City, Country',
                'centre': 'centre1',
                'email': 'e1@mail.com',
                'firstName': 'Mary',
                'laboratory': '',
                'lastName': 'Smith',
                'telephone': ''
            },
            {
                'address': '123 Road Name, City, Country',
                'centre': 'centre1',
                'email': 'e1@mail.com',
                'firstName': 'Robert',
                'laboratory': '',
                'lastName': 'Williams',
                'telephone': ''
            },
            {
                'address': '',
                'centre': 'centre1',
                'email': 'e1@mail.com',
                'firstName': 'Patricia',
                'laboratory': '',
                'lastName': 'Smith',
                'telephone': ''
            },
            {
                'address': '',
                'centre': 'centre1',
                'email': 'e1@mail.com',
                'firstName': 'John',
                'laboratory': '',
                'lastName': 'Brown',
                'telephone': ''
            }
        ]

        self.assertEqual(expected, result)
        mock_fetch.assert_called_once_with("STUDY123")

    # convert_gvf_to_vcf.metadataJSON.DGVaMetadataRetriever.load_from_db                TARGET
    # PATCH                                                                             THE PATCH
    # mock_load                                                                         THE MOCK
    @patch("convert_gvf_to_vcf.metadataJSON.DGVaMetadataRetriever.load_from_db")
    def test__determine_project_pre_registered(self, mock_load):
        # create mock data from load_from_db (PROJECT NOT PRE-REGISTERED)
        mock_data = {
            0: (None,) # not pre-registered
        }
        mock_load.return_value = mock_data

        metadata_client = DGVaMetadataRetriever(self.config)
        is_preregistered, accession = metadata_client._determine_project_pre_registered("STUDY123")
        self.assertFalse(is_preregistered)
        self.assertEqual(accession, None)
        # create mock data from load_from_db (PROJECT PRE-REGISTERED)
        mock_data = {
            0: ('PRJNA28889',) # pre-registered
        }
        mock_load.return_value = mock_data
        metadata_client = DGVaMetadataRetriever(self.config)
        is_preregistered, accession =  metadata_client._determine_project_pre_registered("STUDY789")
        self.assertTrue(is_preregistered)
        self.assertEqual(accession, "PRJNA28889")

    @patch("convert_gvf_to_vcf.metadataJSON.DGVaMetadataRetriever._fetch_reference_genome")
    @patch("convert_gvf_to_vcf.metadataJSON.DGVaMetadataRetriever._fetch_experiment_type")
    @patch("convert_gvf_to_vcf.metadataJSON.DGVaMetadataRetriever._fetch_analysis_description")
    # convert_gvf_to_vcf.metadataJSON.DGVaMetadataRetriever._fetch_analysis_alias                TARGET
    # PATCH                                                                             THE PATCH
    # mock_alias                                                                         THE MOCK
    @patch("convert_gvf_to_vcf.metadataJSON.DGVaMetadataRetriever._fetch_analysis_alias")
    def test__get_analysis(self, mock_alias, mock_desc, mock_experiment_type, mock_reference_genome):
        mock_alias.return_value = "MYanalysisALIAS"
        # mock_title.return_value = "mytitle"
        mock_desc.return_value = "mock_desc"
        mock_experiment_type.return_value = "mock_experiment_type"
        mock_reference_genome.return_value = "GCA000000000"
        metadata_client = DGVaMetadataRetriever(self.config)
        result = metadata_client._get_analysis("estd123", self.vcf_output, "assemble.fasta", "assembly_report.txt")
        expected_result = [{'analysisTitle': '', 'analysisAlias': 'MYanalysisALIAS', 'description': 'mock_desc', 'experimentType': 'mock_experiment_type', 'referenceGenome': 'GCA000000000', 'evidenceType': 'genotype', 'referenceFasta': 'assemble.fasta', 'assemblyReport': 'assembly_report.txt', 'platform': '', 'software': '', 'pipelineDescriptions': '', 'imputation': False, 'phasing': False, 'date': '', 'centre': '', 'links': '', 'runAccessions': []}]
        self.assertEqual(result,expected_result)

    # convert_gvf_to_vcf.metadataJSON.DGVaMetadataRetriever.load_from_db                TARGET
    # PATCH                                                                             THE PATCH
    # mock_load                                                                         THE MOCK
    @patch("convert_gvf_to_vcf.metadataJSON.DGVaMetadataRetriever.load_from_db")
    def test__determine_sample_pre_registered(self, mock_load):
        # (SAMPLE NOT PRE-REGISTERED)
        # create mock data from load_from_db (SAMPLE NOT PRE-REGISTERED)
        mock_data = {
            0: (None, "NA20344"),
            1: (None, "NA20345")  # not pre-registered
        }
        mock_load.return_value = mock_data

        metadata_client = DGVaMetadataRetriever(self.config)
        result = metadata_client._determine_sample_pre_registered("estd123")

        SampleStatus = namedtuple('SampleStatus', ['is_sample_preregistered', 'sample_accession', 'sample_id'])
        expected_result = [SampleStatus(is_sample_preregistered=False, sample_accession=None, sample_id="NA20344"), SampleStatus(is_sample_preregistered=False, sample_accession=None, sample_id="NA20345")]
        self.assertEqual(result, expected_result)

        # create mock data from load_from_db (SAMPLE PRE-REGISTERED)
        mock_data = {
            0: ("SAMN12345678","mypreregisteredsample"),
        }
        mock_load.return_value = mock_data

        metadata_client = DGVaMetadataRetriever(self.config)
        result = metadata_client._determine_sample_pre_registered("estd123")

        SampleStatus = namedtuple('SampleStatus', ['is_sample_preregistered', 'sample_accession', 'sample_id'])
        expected_result = [SampleStatus(is_sample_preregistered=True, sample_accession="SAMN12345678", sample_id="mypreregisteredsample")]
        self.assertEqual(result, expected_result)