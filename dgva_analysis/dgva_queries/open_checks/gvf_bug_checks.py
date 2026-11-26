"""Tally the header keys and column 9 attributes of GVF files.

Runs locally on downloaded GVF files (.gvf or .gvf.gz); it does not connect to the database.
For each file it writes one row per (section, key) with how often the key appears:
  - section 'header:<pragma>' for #subject / #sample / other '#' lines, key = attribute name in that line;
  - section 'call' / 'region' for data lines, key = column 9 attribute name.
It adds derived checks as extra keys:
  - 'Variant_seq:informative'  Variant_seq is not '.', '-' or empty;
  - 'mate_id:equals_own_name'  mate_id equals the line's own Name;
  - 'ciend:left_equals_right'  ciend is written as -X,X.

    python gvf_bug_checks.py --output gvf_bug_checks.tsv file1.gvf.gz ...
"""
import argparse
import csv
import gzip
import os
import re
import sys
from collections import Counter

ACCESSION_KIND = re.compile(r'^(?:e|n)(s)?sv\d+$')


def open_text(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def check_file(path):
    counts = Counter()
    lines = Counter()
    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('#'):
                if line.startswith('##'):
                    continue
                pragma, _, rest = line[1:].partition(':')
                section = f'header:{pragma.strip()}'
                lines[section] += 1
                for item in rest.split(';'):
                    if '=' in item:
                        counts[(section, item.split('=', 1)[0].strip())] += 1
                continue
            columns = line.split('\t')
            attributes = dict(item.split('=', 1) for item in columns[8].split(';') if '=' in item)
            match = ACCESSION_KIND.match(attributes.get('Name', ''))
            section = ('call' if match.group(1) else 'region') if match else 'unknown'
            lines[section] += 1
            for key in attributes:
                counts[(section, key)] += 1
            if attributes.get('Variant_seq', '.') not in ('.', '-', ''):
                counts[(section, 'Variant_seq:informative')] += 1
            if 'mate_id' in attributes and attributes['mate_id'] == attributes.get('Name'):
                counts[(section, 'mate_id:equals_own_name')] += 1
            ciend = attributes.get('ciend')
            if ciend and ciend.lstrip('-').split(',')[0] == ciend.split(',')[-1]:
                counts[(section, 'ciend:left_equals_right')] += 1
    for section, n in lines.items():
        yield section, '(lines)', n
    for (section, key), n in sorted(counts.items()):
        yield section, key, n


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('gvf', nargs='+')
    parser.add_argument('--output', required=True)
    args = parser.parse_args(argv)
    with open(args.output, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t', lineterminator='\n')
        writer.writerow(['file', 'section', 'key', 'count'])
        for path in args.gvf:
            for section, key, n in check_file(path):
                writer.writerow([os.path.basename(path), section, key, n])
            print(f'{path}: done', file=sys.stderr)


if __name__ == '__main__':
    main()
