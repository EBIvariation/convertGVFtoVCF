"""Count element/attribute occurrences in each study's latest-dated Exchange XML in loading_files/.

Reads files only (no database). Usage: python -I scan_exchange.py <paths.tsv> <submissions_dir> <out_dir> [<loaded_files.tsv>]
With loaded_files.tsv (columns study_accession, loaded_file_name, path_in_copy), scans the file at path_in_copy
for each study instead, and writes loaded_fields.tsv.
Writes <out_dir>/exchange_fields.tsv: study, file, path (element or element@attribute), occurrences, non-empty values.
"""
import collections
import csv
import gzip
import os
import re
import sys
import xml.etree.ElementTree as ET

paths_tsv, sub_dir, out_dir = sys.argv[1:4]
loaded = None
if len(sys.argv) > 4:
    loaded = {r['study_accession']: r['path_in_copy'] for r in csv.DictReader(open(sys.argv[4]), delimiter='\t')}
os.makedirs(out_dir, exist_ok=True)


def latest_exchange(folder):
    lf = os.path.join(sub_dir, folder, 'loading_files')
    if not os.path.isdir(lf):
        return None
    cands = []
    for f in os.listdir(lf):
        if re.search(r'exchange\.xml(\.gz)?$', f):
            m = re.search(r'\.(\d{8})\.exchange', f)
            cands.append((m.group(1) if m else '00000000', f))
    if not cands:
        return None
    return os.path.join(lf, max(cands)[1])


rows = list(csv.DictReader(open(paths_tsv), delimiter='\t'))
done = set()
out_path = os.path.join(out_dir, 'loaded_fields.tsv' if loaded else 'exchange_fields.tsv')
if os.path.exists(out_path):
    done = {l.split('\t')[0] for l in open(out_path)}
out = open(out_path, 'a')
for r in rows:
    acc, fid = r['study_accession'], r['folder_id']
    if acc in done:
        continue
    if loaded is not None:
        f = loaded.get(acc) or None
    elif not fid:
        continue
    else:
        f = latest_exchange(fid)
    if not f:
        continue
    counts = collections.Counter()
    nonempty = collections.Counter()
    opener = gzip.open if f.endswith('.gz') else open
    stack = []
    try:
        with opener(f, 'rb') as fh:
            for ev, el in ET.iterparse(fh, events=('start', 'end')):
                if ev == 'start':
                    stack.append(el.tag)
                    continue
                tag = el.tag
                parent = stack[-2] if len(stack) > 1 else ''
                key = f'{parent}/{tag}'
                counts[key] += 1
                if el.text and el.text.strip():
                    nonempty[key] += 1
                for a, v in el.attrib.items():
                    counts[f'{key}@{a}'] += 1
                    if v.strip():
                        nonempty[f'{key}@{a}'] += 1
                stack.pop()
                el.clear()
    except Exception as e:  # keep going; record the error
        out.write(f'{acc}\t{os.path.relpath(f, os.path.dirname(sub_dir))}\tERROR\t0\t{type(e).__name__}: {e}\n')
        out.flush()
        continue
    for k in sorted(counts):
        out.write(f'{acc}\t{os.path.relpath(f, os.path.dirname(sub_dir))}\t{k}\t{counts[k]}\t{nonempty[k]}\n')
    out.flush()
    print(acc, os.path.basename(f), len(counts), flush=True)
