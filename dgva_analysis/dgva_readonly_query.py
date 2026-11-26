"""Run one saved query against the DGVa database through ReadOnlyOracleClient and write the result as TSV.

Example:
    python -m dgva_analysis.dgva_readonly_query --config config.yaml --sql dgva_analysis/dgva_queries/schema/tables.sql \
        --bind owner=DGVA --max-rows 50
"""
import argparse
import csv
import sys

from ebi_eva_common_pyutils.logger import logging_config as log_cfg

from dgva_analysis.read_only_query import ReadOnlyOracleClient, DEFAULT_MAX_ROWS

logger = log_cfg.get_logger(__name__)


def parse_binds(bind_args):
    binds = {}
    for bind in bind_args or []:
        name, sep, value = bind.partition('=')
        if not sep or not name:
            raise argparse.ArgumentTypeError(f'Bind "{bind}" must look like name=value.')
        binds[name] = value
    return binds


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', required=True, help='YAML file with a DGVA section (host, port, user, password, '
                                                         'service_name).')
    parser.add_argument('--sql', required=True, help='File holding one SELECT or WITH statement.')
    parser.add_argument('--bind', action='append', help='Bind variable as name=value (repeatable).')
    parser.add_argument('--max-rows', type=int, default=DEFAULT_MAX_ROWS,
                        help=f'Row limit (default {DEFAULT_MAX_ROWS}). 0 means no limit.')
    parser.add_argument('--output', help='TSV output file (default: standard output).')
    args = parser.parse_args(argv)
    log_cfg.add_stdout_handler()

    max_rows = None if args.max_rows == 0 else args.max_rows
    with ReadOnlyOracleClient.from_config(args.config) as client:
        result = client.query_file(args.sql, parse_binds(args.bind), max_rows=max_rows)

    out = open(args.output, 'w', newline='') if args.output else sys.stdout
    try:
        writer = csv.writer(out, delimiter='\t', lineterminator='\n')
        writer.writerow(result.columns)
        writer.writerows(result.rows)
    finally:
        if args.output:
            out.close()
    logger.info(f'{len(result.rows)} rows written' + (' (cut by the row limit)' if result.truncated else ''))


if __name__ == '__main__':
    main()
