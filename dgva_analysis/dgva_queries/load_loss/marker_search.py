"""Count, in every Exchange XML under <dgva_root>/submissions and <dgva_root>/studies, the GENOTYPE and ALLELE elements
and the non-empty mode_of_inheritance, inheritance, event ID, mate ID and variant_region_group_id attributes.
Reads files only.

Usage: python -I marker_search.py <dgva_root> <out.tsv>
"""
import collections
import gzip
import multiprocessing
import os
import re
import sys

PAT = re.compile(rb'<GENOTYPE[ >]|<ALLELE[ >]|\b(mode_of_inheritance|inheritance|variant_call_event_id|'
                 rb'variant_region_event_id|variant_call_mate_id|variant_region_mate_id|variant_region_group_id|'
                 rb'event_id|mate_id)="[^"]+"')


def scan(path):
    c = collections.Counter()
    try:
        with open(path, 'rb') as fh:
            magic = fh.read(2)
        fh = gzip.open(path, 'rb') if magic == b'\x1f\x8b' else open(path, 'rb')
        with fh:
            for line in fh:
                for m in PAT.finditer(line):
                    c[(m.group(1) or m.group(0)[1:-1]).decode()] += 1
    except Exception as e:
        c['ERROR:' + type(e).__name__] += 1
    return path, c


if __name__ == '__main__':
    root, out = sys.argv[1:3]
    files = []
    for top in ('submissions', 'studies'):
        for r, ds, fs in os.walk(os.path.join(root, top)):
            for f in fs:
                if 'exchange.xml' in f and 'ctl' not in f and 'md5' not in f:
                    files.append(os.path.join(r, f))
    with multiprocessing.Pool(6) as pool, open(out, 'w') as o:
        for path, c in pool.imap_unordered(scan, files, chunksize=4):
            if c:
                o.write(os.path.relpath(path, root) + '\t' + ';'.join(f'{k}={v}' for k, v in sorted(c.items())) + '\n')
                o.flush()
    print('files scanned:', len(files))
