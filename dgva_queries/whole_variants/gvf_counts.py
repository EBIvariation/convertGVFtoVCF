"""Ticket 6: count the variants in FTP GVF files, to compare with exporter_drops.sql.

Runs locally on downloaded GVF files (.gvf or .gvf.gz); it does not connect to the database.
For each file it counts, separately for calls (essv/nssv) and regions (esv/nsv):
  - lines, and distinct accessions (the Name attribute); one call can be written on several lines (bug B1);
  - lines with no Name.
The study accession, assembly name and label are read from the file name:
  <accession>_<display name>.<date>.<assembly>.<Submitted|Remapped>.gvf[.gz]

    python gvf_counts.py --output gvf_counts.tsv file1.gvf.gz file2.gvf.gz ...

With --names-dir, it also writes the distinct accessions of each file to <names-dir>/<file>.names, one per line.
"""
import argparse
import csv
import gzip
import os
import re
import sys
from collections import defaultdict

FILE_NAME = re.compile(r'^(?P<study>[en]std\d+)_.*\.(?P<date>\d{4}-\d{2}-\d{2})\.(?P<assembly>.+)\.'
                       r'(?P<label>Submitted|Remapped)\.gvf(?:\.gz)?$')
ACCESSION_KIND = re.compile(r'^(?:e|n)(?P<call>s)?sv\d+$')


def open_text(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def count_file(path, names_dir=None):
    match = FILE_NAME.match(os.path.basename(path))
    if not match:
        raise ValueError(f'Unexpected file name: {path}')
    lines = defaultdict(int)
    names = defaultdict(set)
    with open_text(path) as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            attributes = dict(item.split('=', 1) for item in line.rstrip('\n').split('\t')[8].split(';') if '=' in item)
            name = attributes.get('Name', '')
            kind_match = ACCESSION_KIND.match(name)
            kind = ('call' if kind_match.group('call') else 'region') if kind_match else 'no_or_unknown_name'
            lines[kind] += 1
            if kind_match:
                names[kind].add(name)
    if names_dir:
        os.makedirs(names_dir, exist_ok=True)
        with open(os.path.join(names_dir, os.path.basename(path) + '.names'), 'w') as out:
            for name in sorted(names['call'] | names['region']):
                out.write(name + '\n')
    for kind in sorted(lines):
        yield {**match.groupdict(), 'file': os.path.basename(path), 'region_or_call': kind,
               'lines': lines[kind], 'distinct_accessions': len(names[kind])}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('gvf', nargs='+')
    parser.add_argument('--output', required=True)
    parser.add_argument('--names-dir')
    args = parser.parse_args(argv)
    columns = ['study', 'assembly', 'label', 'date', 'region_or_call', 'lines', 'distinct_accessions', 'file']
    with open(args.output, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=columns, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for path in args.gvf:
            for row in count_file(path, args.names_dir):
                writer.writerow(row)
            print(f'{path}: done', file=sys.stderr)


if __name__ == '__main__':
    main()
