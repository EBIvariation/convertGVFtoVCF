from abc import ABC, abstractmethod

import oracledb
from ebi_eva_common_pyutils.config import cfg
from ebi_eva_common_pyutils.logger import logging_config as log_cfg

logger = log_cfg.get_logger(__name__)
class BaseMetadataRetriever(ABC):
    """ The responsibility of this base class is to ensure a consistent method for metadata retrieval.
    This includes database connection and the ingestion and validation of data.
    """
    def __init__(self, path_to_config_yaml):
        # coming from the config file
        cfg.load_config_file(path_to_config_yaml)  # cfg is a dictionary
        # db connection setup
        self._connection = None
        self._host = self._get_validated_value(cfg, ("DGVA","host"), str, default_value=None) # get information from the config dictionary
        self._port = self._get_validated_value(cfg, ("DGVA", "port"), str, default_value=None)
        self._username = self._get_validated_value(cfg, ("DGVA", "user"), str, default_value=None)
        self._password = self._get_validated_value(cfg, ("DGVA", "password"), str, default_value=None)
        self._service_name = self._get_validated_value(cfg, ("DGVA", "service_name"), str, default_value=None)
        # db parameters
        self._max_retries = 3

    def __enter__(self):
        # enables the "with" context manager
        _ = self.connection
        logger.info("Opening connection using a context manager - SUCCESS")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # closes the connection and resets the state of connection
        if self._connection:
            self._connection.close()
            self._connection = None
            logger.info("Closing connection safely using a context manager - SUCCESS")

    @staticmethod
    def _get_validated_value(config, key_parts_to_get, expected_type, default_value=None):
        # key_parts_to_get is a tuple of multiple values so unpacking
        value_in_config = config.query(*key_parts_to_get, ret_default=default_value)
        if value_in_config is None:
            raise ValueError(f"Missing key required: {key_parts_to_get}")
        if not isinstance(value_in_config, expected_type):
            try:
                # cast the value to its expected type
                value_in_config = expected_type(value_in_config)
            except (ValueError, TypeError):
                raise TypeError(f"Key '{key_parts_to_get}' must be {expected_type.__name__}")
        return value_in_config

    @property
    def connection(self):
        attempts = 0
        while attempts < self._max_retries:
            # if there is a healthy connection established, return it
            if self._connection is not None and self._connection.is_healthy():
                logger.info("A healthy connection has already been established.")
                return self._connection
            # if there is an unhealthy connection, close it and reset the state
            if self._connection:
                logger.warning("Unhealthy connection. Attempting to reconnect.")
                try:
                    # close the unhealthy connection
                    self._connection.close()
                    logger.info("Connection closed safely.")
                except:
                    pass
                self._connection = None
            # first connection or reconnection
            attempts += 1
            try:
                logger.info("Connecting to the database.")
                self._connection = oracledb.connect(host=self._host, port=self._port,
                                                    user=self._username, password=self._password, service_name=self._service_name)
                return self._connection
            except Exception as e:
                logger.error(f"Failed to connect. {e}")
                if attempts >= self._max_retries:
                    logger.error(f"Max connection retries reached: {self._max_retries}.")
                    raise

    def load_from_db(self, query_to_load):
        """ Searches the DGVa database using the query provided.
        :param query_to_load : SQL query
        :return: a list of rows (or if query returns nothing an empty list)
        """
        try:
            # create the iterator to process queries
            with self.connection.cursor() as cur:
                # prepare to fetch data
                ###############################################################################
                # WARNING: do not use f-strings or % formatting here, use placeholders instead
                query = query_to_load
                # query_placeholders = query_place_holder_to_load
                # cur.execute(query, query_placeholders)
                ###############################################################################
                # run the sql query
                logger.info(f"Executing the following query: {query}")
                cur.execute(query)
                # fetch the data from the cursor
                rows = cur.fetchall() # returns list of tuples
                logger.info(f"Fetching the data: {rows}")
                logger.info(f"Fetching metadata query - SUCCESS - {len(rows)} records found")
                return rows
        except Exception as e:
            logger.warning(f"Database error: {e}")
            return []

    def validate_fetch_result(self, eva_field_name, fetch_result_list, single_value_as_str):
        try:
            if fetch_result_list:
                results = [row[0] for row in fetch_result_list if row]
                if single_value_as_str:
                    fetch_result = results[0] if len(results) == 1 else results
                else:
                    fetch_result = results
                if fetch_result:
                    # SUCCESS if value is present or None
                    logger.info(f"Fetching {eva_field_name} - SUCCESS - Value(s) for {eva_field_name} found: {fetch_result}.")
            else:
                raise ValueError(f"Missing data: {eva_field_name}.")
        except ValueError as e:
            logger.error(f"Fetching {eva_field_name}  - FAILURE - {eva_field_name} not found. {e} Setting value as empty string.")
            fetch_result = ""
        return fetch_result

    @abstractmethod
    def retrieve(self):
        pass