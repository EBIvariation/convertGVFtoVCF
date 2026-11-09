"""Read-only access to the DGVa Oracle database.

The DGVa account can write, and the database is the only live copy of DGVa. Every query must therefore go through
ReadOnlyOracleClient, which:

- checks each statement before anything is sent to the database (check_read_only): one statement, starting with
  SELECT or WITH, with no DML, DDL, PL/SQL, transaction control, SELECT ... FOR UPDATE, database link, or call to a
  function that is not a known Oracle built-in;
- runs each statement in its own transaction opened with SET TRANSACTION READ ONLY;
- never commits, and rolls back after every statement and before closing;
- refuses a connection with autocommit on;
- takes values as bind variables, and applies a row limit unless told otherwise.

The checks are a second line of defence. An account with SELECT grants only is still the only safeguard that also
blocks DDL.
"""
import re
from dataclasses import dataclass

import oracledb
from ebi_eva_common_pyutils.config import cfg
from ebi_eva_common_pyutils.logger import logging_config as log_cfg

from convert_gvf_to_vcf.utils import get_validated_value

logger = log_cfg.get_logger(__name__)

DEFAULT_MAX_ROWS = 1000
READ_ONLY_TRANSACTION = 'SET TRANSACTION READ ONLY'

# Words that never appear in a read-only query. Matched as whole words, outside string literals, quoted identifiers
# and comments, so column names such as UPDATE_COMMENT or STUDY_FILE.FILE_TYPE are not affected.
FORBIDDEN_WORDS = {
    # DML
    'INSERT', 'UPDATE', 'DELETE', 'MERGE', 'UPSERT', 'INTO',
    # DDL and privileges
    'CREATE', 'ALTER', 'DROP', 'TRUNCATE', 'RENAME', 'GRANT', 'REVOKE', 'COMMENT', 'PURGE', 'FLASHBACK', 'ANALYZE',
    'AUDIT', 'NOAUDIT', 'ASSOCIATE', 'DISASSOCIATE',
    # transaction and session control
    'COMMIT', 'ROLLBACK', 'SAVEPOINT', 'LOCK', 'SET',
    # PL/SQL
    'BEGIN', 'DECLARE', 'EXEC', 'EXECUTE', 'CALL', 'FUNCTION', 'PROCEDURE', 'PACKAGE', 'TRIGGER',
}

# SQL keywords that can be followed by an opening parenthesis without being a function call.
KEYWORDS_BEFORE_PARENTHESIS = {
    'SELECT', 'FROM', 'WHERE', 'AND', 'OR', 'NOT', 'IN', 'EXISTS', 'ANY', 'ALL', 'SOME', 'AS', 'ON', 'USING', 'JOIN',
    'BY', 'GROUP', 'HAVING', 'UNION', 'INTERSECT', 'MINUS', 'EXCEPT', 'CASE', 'WHEN', 'THEN', 'ELSE', 'BETWEEN', 'LIKE',
    'DISTINCT', 'OVER', 'WITHIN', 'KEEP', 'PIVOT', 'UNPIVOT', 'FOR', 'LATERAL', 'ESCAPE', 'PARTITION', 'IS', 'WITH',
}

# Oracle built-in functions that only read. Anything else followed by "(" is rejected, because a user-defined
# function can write through an autonomous transaction even inside a read-only transaction.
# Add a function here only after checking it cannot write.
ALLOWED_FUNCTIONS = {
    # aggregates and analytics
    'COUNT', 'SUM', 'MIN', 'MAX', 'AVG', 'MEDIAN', 'STDDEV', 'VARIANCE', 'LISTAGG', 'ROW_NUMBER', 'RANK',
    'DENSE_RANK', 'LAG', 'LEAD', 'FIRST_VALUE', 'LAST_VALUE', 'NTILE', 'PERCENTILE_CONT', 'PERCENTILE_DISC',
    'APPROX_COUNT_DISTINCT',
    # null handling and comparison
    'COALESCE', 'NVL', 'NVL2', 'NULLIF', 'DECODE', 'GREATEST', 'LEAST', 'LNNVL',
    # conversion
    'CAST', 'TO_CHAR', 'TO_NUMBER', 'TO_DATE', 'TO_TIMESTAMP', 'TO_CLOB', 'EXTRACT', 'NUMTODSINTERVAL',
    'NUMTOYMINTERVAL',
    # strings
    'UPPER', 'LOWER', 'INITCAP', 'TRIM', 'LTRIM', 'RTRIM', 'LPAD', 'RPAD', 'SUBSTR', 'INSTR', 'LENGTH', 'REPLACE',
    'TRANSLATE', 'CONCAT', 'REGEXP_LIKE', 'REGEXP_SUBSTR', 'REGEXP_INSTR', 'REGEXP_REPLACE', 'REGEXP_COUNT', 'ASCII',
    'CHR', 'STANDARD_HASH', 'ORA_HASH',
    # numbers and dates
    'ABS', 'ROUND', 'TRUNC', 'FLOOR', 'CEIL', 'MOD', 'SIGN', 'POWER', 'SQRT', 'LOG', 'LN', 'EXP', 'SYSDATE',
    'ADD_MONTHS', 'MONTHS_BETWEEN', 'LAST_DAY',
}
# Package functions that only read (for CLOB columns such as VARIANT_CALL_SEQUENCE.VARIANT_CALL_SEQ).
ALLOWED_QUALIFIED_FUNCTIONS = {'DBMS_LOB.GETLENGTH', 'DBMS_LOB.SUBSTR', 'DBMS_LOB.INSTR'}

_Q_QUOTE_CLOSERS = {'[': ']', '{': '}', '(': ')', '<': '>'}


class UnsafeQueryError(ValueError):
    """Raised when a statement is not provably read-only. Nothing has been sent to the database."""


@dataclass
class QueryResult:
    columns: list
    rows: list
    truncated: bool  # True when the row limit cut the result short


def _mask(sql):
    """Return the statement with comments removed, string literals replaced by '' and quoted identifiers replaced by
    the word QUOTED_IDENTIFIER, so that keywords can be searched for safely."""
    out = []
    i, n = 0, len(sql)
    while i < n:
        ch = sql[i]
        nxt = sql[i + 1] if i + 1 < n else ''
        if ch == '-' and nxt == '-':
            end = sql.find('\n', i)
            i = n if end == -1 else end
            out.append(' ')
        elif ch == '/' and nxt == '*':
            end = sql.find('*/', i + 2)
            if end == -1:
                raise UnsafeQueryError('Unterminated comment.')
            i = end + 2
            out.append(' ')
        elif ch in 'qQ' and nxt == "'" and (i == 0 or not (sql[i - 1].isalnum() or sql[i - 1] == '_')):
            # Oracle alternative quoting: q'[...]', q'{...}', q'X...X'
            if i + 2 >= n:
                raise UnsafeQueryError('Unterminated string literal.')
            opener = sql[i + 2]
            closer = _Q_QUOTE_CLOSERS.get(opener, opener)
            end = sql.find(closer + "'", i + 3)
            if end == -1:
                raise UnsafeQueryError('Unterminated string literal.')
            i = end + 2
            out.append("''")
        elif ch == "'":
            j = i + 1
            while True:
                end = sql.find("'", j)
                if end == -1:
                    raise UnsafeQueryError('Unterminated string literal.')
                if end + 1 < n and sql[end + 1] == "'":  # escaped quote
                    j = end + 2
                    continue
                break
            i = end + 1
            out.append("''")
        elif ch == '"':
            end = sql.find('"', i + 1)
            if end == -1:
                raise UnsafeQueryError('Unterminated quoted identifier.')
            i = end + 1
            out.append(' QUOTED_IDENTIFIER ')
        else:
            out.append(ch)
            i += 1
    return ''.join(out)


def check_read_only(sql):
    """Check that a statement is a single read-only query. Raise UnsafeQueryError otherwise.

    :param sql: the statement, optionally ending with one semicolon
    :return: the statement to execute (trailing semicolon and surrounding whitespace removed)
    """
    if not isinstance(sql, str) or not sql.strip():
        raise UnsafeQueryError('Empty statement.')
    statement = sql.strip()
    masked = _mask(statement).strip()
    if masked.endswith(';'):
        masked = masked[:-1].rstrip()
        statement = statement.rstrip()[:-1].rstrip()
    if ';' in masked:
        raise UnsafeQueryError('Only one statement is allowed.')
    if '@' in masked:
        raise UnsafeQueryError('Database links (@) are not allowed.')

    upper = masked.upper()
    first_word = re.match(r'\s*([A-Z_]+)', upper)
    if not first_word or first_word.group(1) not in ('SELECT', 'WITH'):
        raise UnsafeQueryError('The statement must start with SELECT or WITH.')
    if re.search(r'\bFOR\s+UPDATE\b', upper):
        raise UnsafeQueryError('SELECT ... FOR UPDATE is not allowed.')
    words = set(re.findall(r'\b[A-Z_][A-Z0-9_$#]*\b', upper))
    forbidden = sorted(words & FORBIDDEN_WORDS)
    if forbidden:
        raise UnsafeQueryError(f'Forbidden keyword(s): {", ".join(forbidden)}.')

    # Every "name(" must be a keyword, an allowed built-in function, or a WITH-clause name.
    cte_names = set(re.findall(r'(?:\bWITH|,)\s*([A-Z_][A-Z0-9_$#]*)\s*(?:\([^()]*\)\s*)?AS\s*\(', upper))
    for match in re.finditer(r'([A-Z_][A-Z0-9_$#]*(?:\s*\.\s*[A-Z_][A-Z0-9_$#]*)*)\s*\(', upper):
        name = re.sub(r'\s+', '', match.group(1))
        if '.' in name:
            if name not in ALLOWED_QUALIFIED_FUNCTIONS:
                raise UnsafeQueryError(f'Call to {name} is not allowed (PL/SQL package or schema function).')
        elif name not in KEYWORDS_BEFORE_PARENTHESIS | ALLOWED_FUNCTIONS | cte_names:
            raise UnsafeQueryError(f'Call to {name} is not allowed: it is not a known read-only built-in function.')
    return statement


class ReadOnlyOracleClient:
    """Runs read-only queries against the DGVa Oracle database. Use it as a context manager:

        with ReadOnlyOracleClient.from_config('config.yaml') as client:
            result = client.query('SELECT display_name FROM dgva.dgva_study WHERE study_accession = :acc',
                                  {'acc': 'estd1'}, max_rows=10)
    """

    def __init__(self, host, port, user, password, service_name, module='convertGVFtoVCF-readonly'):
        self._connect_args = dict(host=host, port=port, user=user, password=password, service_name=service_name)
        self._module = module
        self._connection = None

    @classmethod
    def from_config(cls, path_to_config_yaml):
        """Read the connection details from the DGVA section of a YAML config file, as BaseMetadataRetriever does."""
        cfg.load_config_file(path_to_config_yaml)
        return cls(
            host=get_validated_value(cfg, ('DGVA', 'host'), str, default_value=None),
            port=get_validated_value(cfg, ('DGVA', 'port'), str, default_value=None),
            user=get_validated_value(cfg, ('DGVA', 'user'), str, default_value=None),
            password=get_validated_value(cfg, ('DGVA', 'password'), str, default_value=None),
            service_name=get_validated_value(cfg, ('DGVA', 'service_name'), str, default_value=None),
        )

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def _get_connection(self):
        if self._connection is None:
            logger.info(f'Connecting read-only to {self._connect_args["host"]}/{self._connect_args["service_name"]}')
            connection = oracledb.connect(**self._connect_args)
            if connection.autocommit:
                connection.close()
                raise RuntimeError('The connection has autocommit on. Refusing to use it.')
            connection.module = self._module
            self._connection = connection
        return self._connection

    def query(self, sql, binds=None, max_rows=DEFAULT_MAX_ROWS):
        """Run one read-only query.

        :param sql: a SELECT or WITH statement; pass values as bind variables (:name), never by string formatting
        :param binds: dict (or sequence) of bind values
        :param max_rows: maximum number of rows to return; None for no limit
        :return: QueryResult
        """
        statement = check_read_only(sql)  # raises before anything is sent to the database
        if binds is None:
            binds = {}
        binds = dict(binds) if isinstance(binds, dict) else list(binds)
        if max_rows is not None:
            if not isinstance(max_rows, int) or max_rows < 1:
                raise ValueError('max_rows must be a positive integer or None.')
            if not isinstance(binds, dict):
                raise ValueError('Use named bind variables when max_rows is set.')
            statement = f'SELECT * FROM (\n{statement}\n) WHERE ROWNUM <= :readonly_helper_max_rows'
            binds['readonly_helper_max_rows'] = max_rows + 1

        connection = self._get_connection()
        try:
            connection.rollback()  # end any open transaction, so that the next statement starts a new one
            with connection.cursor() as cursor:
                cursor.execute(READ_ONLY_TRANSACTION)
                cursor.execute(statement, binds)
                columns = [d[0] for d in cursor.description]
                rows = cursor.fetchall()
        finally:
            connection.rollback()
        truncated = max_rows is not None and len(rows) > max_rows
        if truncated:
            rows = rows[:max_rows]
            logger.warning(f'Result cut at {max_rows} rows.')
        return QueryResult(columns=columns, rows=rows, truncated=truncated)

    def query_file(self, path, binds=None, max_rows=DEFAULT_MAX_ROWS):
        """Run the single query saved in a .sql file."""
        with open(path) as sql_file:
            return self.query(sql_file.read(), binds, max_rows)

    def close(self):
        if self._connection is not None:
            try:
                self._connection.rollback()
            finally:
                self._connection.close()
                self._connection = None
                logger.info('Read-only connection closed.')
