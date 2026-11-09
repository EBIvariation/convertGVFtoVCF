import os
from unittest import TestCase
from unittest.mock import patch, MagicMock, call

from convert_gvf_to_vcf.metadata_retrievers.read_only_query import (
    check_read_only, ReadOnlyOracleClient, UnsafeQueryError, READ_ONLY_TRANSACTION)

CONNECT = 'convert_gvf_to_vcf.metadata_retrievers.read_only_query.oracledb.connect'

# Statements that must be rejected, grouped by kind. None of them may reach the database.
UNSAFE = {
    'insert': [
        "INSERT INTO dgva.study_comment (study_comment) VALUES ('x')",
        "insert all into t values (1) select * from dual",
    ],
    'update': [
        "UPDATE dgva.dgva_study SET display_name = 'x' WHERE study_accession = 'estd1'",
        "  update dgva_study set hold_date = null",
    ],
    'delete': [
        "DELETE FROM dgva.variant_withdrawn",
        "DELETE dgva.variant_withdrawn WHERE 1 = 1",
    ],
    'merge': [
        "MERGE INTO dgva.dgva_study d USING dual ON (1 = 1) WHEN MATCHED THEN UPDATE SET d.display_name = 'x'",
    ],
    'ddl': [
        "CREATE TABLE estd1_vcf_file_vc AS SELECT * FROM dgva.variant_call",
        "CREATE GLOBAL TEMPORARY TABLE tmp (a NUMBER)",
        "CREATE INDEX i ON dgva.variant_call (variant_call_accession)",
        "CREATE OR REPLACE VIEW v AS SELECT 1 x FROM dual",
        "DROP TABLE estd1_vcf_file_vc",
        "ALTER TABLE dgva.dgva_study ADD (x NUMBER)",
        "TRUNCATE TABLE dgva.study_comment",
        "RENAME a TO b",
        "GRANT SELECT ON dgva.dgva_study TO public",
        "REVOKE SELECT ON dgva.dgva_study FROM public",
        "COMMENT ON TABLE dgva.dgva_study IS 'x'",
        "PURGE RECYCLEBIN",
        "ANALYZE TABLE dgva.dgva_study COMPUTE STATISTICS",
    ],
    'plsql': [
        "BEGIN DBMS_STATS.GATHER_TABLE_STATS('dgva', 'ESTD1_VCF_FILE_VC'); END;",
        "DECLARE x NUMBER; BEGIN NULL; END;",
        "CALL dgva.delete_partof_study_pkg.delete_study('estd1')",
        "EXEC dgva.accession_generator_pkg.next_accession",
        "EXECUTE IMMEDIATE 'DROP TABLE x'",
        "WITH FUNCTION f RETURN NUMBER IS BEGIN RETURN 1; END; SELECT f FROM dual",
        "SELECT dgva.accession_generator_pkg.next_accession() FROM dual",
        "SELECT DBMS_STATS.GATHER_TABLE_STATS('dgva', 'X') FROM dual",
        "SELECT my_function(study_accession) FROM dgva.dgva_study",
        'SELECT "MY_FUNCTION"(study_accession) FROM dgva.dgva_study',
        "SELECT * FROM TABLE(dgva.some_pipelined_function())",
    ],
    'for_update': [
        "SELECT * FROM dgva.dgva_study FOR UPDATE",
        "SELECT * FROM dgva.dgva_study WHERE study_accession = :a FOR UPDATE NOWAIT",
        "select * from dgva.dgva_study for\n  update of display_name",
    ],
    'transaction_control': [
        "COMMIT",
        "ROLLBACK",
        "SAVEPOINT a",
        "SET TRANSACTION READ WRITE",
        "LOCK TABLE dgva.dgva_study IN EXCLUSIVE MODE",
        "ALTER SESSION SET current_schema = dgva",
    ],
    'several_statements': [
        "SELECT 1 FROM dual; DELETE FROM dgva.study_comment",
        "SELECT 1 FROM dual; SELECT 2 FROM dual",
    ],
    'other': [
        "",
        "   ",
        "-- only a comment",
        "EXPLAIN PLAN FOR SELECT * FROM dual",
        "SELECT * INTO x FROM dual",
        "SELECT * FROM dgva.dgva_study@other_db",
        "SELECT 'unterminated FROM dual",
        "SELECT 1 FROM dual /* unterminated",
        "(SELECT 1 FROM dual)",
    ],
}

SAFE = [
    "SELECT display_name FROM dgva.dgva_study WHERE study_accession = :acc",
    "select count(*) from dgva.variant_study where study_accession = :acc and region_or_call = 'call'",
    "SELECT 1 FROM dual;",
    "SELECT 1 FROM dual ;  \n",
    # forbidden words inside literals, quoted identifiers and comments are fine
    "SELECT 'DELETE FROM x; DROP TABLE y' AS s FROM dual",
    "SELECT q'[it's an UPDATE; really]' FROM dual",
    "SELECT \"DROP\" FROM dual",
    "-- DROP TABLE x\nSELECT 1 FROM dual /* COMMIT */",
    "SELECT 'a@b.org' AS email FROM dual",
    # column names containing forbidden words
    "SELECT update_comment, new_feature, correction FROM dgva.study_update WHERE study_accession = :acc",
    "SELECT study_comment, comment_timestamp FROM dgva.study_comment",
    "SELECT file_type, date_added FROM dgva.study_file",
    # joins, CTEs, built-in functions and analytics
    "SELECT * FROM dgva.variant_placement JOIN dgva.variant_study USING (variant_study_id) WHERE study_accession = :acc",
    "WITH calls AS (SELECT variant_accession FROM dgva.variant_study WHERE region_or_call = 'call'), "
    "regions (acc) AS (SELECT variant_accession FROM dgva.variant_study WHERE region_or_call = 'region') "
    "SELECT COUNT(*) FROM calls WHERE variant_accession NOT IN (SELECT acc FROM regions)",
    "SELECT COALESCE(s.start_pos, COALESCE(o.outer_start_pos, i.inner_start_pos)) AS start_pos FROM dual",
    "SELECT ROW_NUMBER() OVER (PARTITION BY study_accession ORDER BY created_date) rn FROM dgva.variant_study",
    "SELECT DBMS_LOB.GETLENGTH(variant_call_seq), DBMS_LOB.SUBSTR(variant_call_seq, 100, 1) "
    "FROM dgva.variant_call_sequence",
    "SELECT * FROM dgva.dgva_study WHERE EXISTS (SELECT 1 FROM dgva.study_remap r WHERE r.study_accession = :acc)",
    "SELECT LISTAGG(taxonomy_id, ',') WITHIN GROUP (ORDER BY taxonomy_id) FROM dgva.study_organism",
    "SELECT * FROM \"DGVA\".\"ALIGNMENT_CV\"",
]


class TestCheckReadOnly(TestCase):

    def test_unsafe_statements_are_rejected(self):
        for kind, statements in UNSAFE.items():
            for sql in statements:
                with self.subTest(kind=kind, sql=sql):
                    with self.assertRaises(UnsafeQueryError):
                        check_read_only(sql)

    def test_safe_statements_are_accepted(self):
        for sql in SAFE:
            with self.subTest(sql=sql):
                check_read_only(sql)

    def test_trailing_semicolon_is_removed(self):
        assert check_read_only("SELECT 1 FROM dual ;  \n") == "SELECT 1 FROM dual"
        assert check_read_only("SELECT ';' FROM dual;") == "SELECT ';' FROM dual"


class TestReadOnlyOracleClient(TestCase):

    def setUp(self):
        self.config = os.path.join(os.path.dirname(__file__), 'input', 'test.config')
        self.connection = MagicMock()
        self.connection.autocommit = False
        self.cursor = self.connection.cursor.return_value.__enter__.return_value
        self.cursor.description = [('STUDY_ACCESSION',), ('DISPLAY_NAME',)]
        self.cursor.fetchall.return_value = [('estd1', 'Redon_et_al_2006')]

    def test_from_config(self):
        client = ReadOnlyOracleClient.from_config(self.config)
        assert client._connect_args == dict(host='mockhost', port='1111', user='mockuser', password='mockpassword',
                                            service_name='mockservicename')

    @patch(CONNECT)
    def test_unsafe_statements_never_reach_the_database(self, mock_connect):
        mock_connect.return_value = self.connection
        with ReadOnlyOracleClient.from_config(self.config) as client:
            for statements in UNSAFE.values():
                for sql in statements:
                    with self.subTest(sql=sql):
                        with self.assertRaises(UnsafeQueryError):
                            client.query(sql)
        # rejected before connecting: no connection was opened and nothing was executed
        mock_connect.assert_not_called()
        self.cursor.execute.assert_not_called()
        self.connection.commit.assert_not_called()

    @patch(CONNECT)
    def test_unsafe_statement_on_open_connection_is_not_executed(self, mock_connect):
        mock_connect.return_value = self.connection
        with ReadOnlyOracleClient.from_config(self.config) as client:
            client.query('SELECT 1 FROM dual')
            self.cursor.execute.reset_mock()
            with self.assertRaises(UnsafeQueryError):
                client.query("SELECT * FROM dgva.dgva_study FOR UPDATE")
            self.cursor.execute.assert_not_called()

    @patch(CONNECT)
    def test_query_runs_in_a_read_only_transaction(self, mock_connect):
        mock_connect.return_value = self.connection
        with ReadOnlyOracleClient.from_config(self.config) as client:
            result = client.query('SELECT study_accession, display_name FROM dgva.dgva_study '
                                  'WHERE study_accession = :acc', {'acc': 'estd1'}, max_rows=None)
        executed = self.cursor.execute.call_args_list
        assert executed[0] == call(READ_ONLY_TRANSACTION)
        assert executed[1] == call('SELECT study_accession, display_name FROM dgva.dgva_study '
                                   'WHERE study_accession = :acc', {'acc': 'estd1'})
        assert len(executed) == 2
        # rolled back before the query, after it, and before closing; never committed
        assert self.connection.rollback.call_count == 3
        self.connection.commit.assert_not_called()
        self.connection.close.assert_called_once()
        assert result.columns == ['STUDY_ACCESSION', 'DISPLAY_NAME']
        assert result.rows == [('estd1', 'Redon_et_al_2006')]
        assert not result.truncated

    @patch(CONNECT)
    def test_rollback_happens_when_the_query_fails(self, mock_connect):
        mock_connect.return_value = self.connection
        self.cursor.execute.side_effect = [None, RuntimeError('ORA-00942')]
        with ReadOnlyOracleClient.from_config(self.config) as client:
            with self.assertRaises(RuntimeError):
                client.query('SELECT * FROM dgva.missing_table')
            assert self.connection.rollback.call_count == 2
        self.connection.commit.assert_not_called()

    @patch(CONNECT)
    def test_row_limit(self, mock_connect):
        mock_connect.return_value = self.connection
        self.cursor.fetchall.return_value = [(i, str(i)) for i in range(3)]
        with ReadOnlyOracleClient.from_config(self.config) as client:
            result = client.query('SELECT a, b FROM t WHERE c = :c;', {'c': 1}, max_rows=2)
        sql, binds = self.cursor.execute.call_args_list[1].args
        assert sql == 'SELECT * FROM (\nSELECT a, b FROM t WHERE c = :c\n) WHERE ROWNUM <= :readonly_helper_max_rows'
        assert binds == {'c': 1, 'readonly_helper_max_rows': 3}
        assert result.rows == [(0, '0'), (1, '1')]
        assert result.truncated

    @patch(CONNECT)
    def test_invalid_row_limit(self, mock_connect):
        mock_connect.return_value = self.connection
        with ReadOnlyOracleClient.from_config(self.config) as client:
            for max_rows in (0, -1, 1.5):
                with self.assertRaises(ValueError):
                    client.query('SELECT 1 FROM dual', max_rows=max_rows)
        mock_connect.assert_not_called()

    @patch(CONNECT)
    def test_autocommit_connection_is_refused(self, mock_connect):
        self.connection.autocommit = True
        mock_connect.return_value = self.connection
        with ReadOnlyOracleClient.from_config(self.config) as client:
            with self.assertRaises(RuntimeError):
                client.query('SELECT 1 FROM dual')
        self.cursor.execute.assert_not_called()
        self.connection.close.assert_called_once()

    @patch(CONNECT)
    def test_query_file(self, mock_connect):
        mock_connect.return_value = self.connection
        sql_file = os.path.join(os.path.dirname(__file__), '..', 'dgva_queries', 'schema', 'tables.sql')
        with ReadOnlyOracleClient.from_config(self.config) as client:
            client.query_file(sql_file, {'owner': 'DGVA'}, max_rows=10)
        sql, binds = self.cursor.execute.call_args_list[1].args
        assert 'FROM all_tables' in sql
        assert binds == {'owner': 'DGVA', 'readonly_helper_max_rows': 11}

    def test_saved_queries_pass_the_check(self):
        query_dir = os.path.join(os.path.dirname(__file__), '..', 'dgva_queries')
        for root, _, files in os.walk(query_dir):
            for name in files:
                if name.endswith('.sql'):
                    with self.subTest(file=name), open(os.path.join(root, name)) as f:
                        check_read_only(f.read())
