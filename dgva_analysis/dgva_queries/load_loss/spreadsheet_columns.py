"""List the columns of DGVa submitter spreadsheets that the spreadsheet-to-Exchange converter never reads (files only).

Usage: PYTHONPATH=<dir with xlrd> python -P -s spreadsheet_columns.py <dgva_root> <dgva-pipeline dir> <study_storage_paths.tsv> <out.tsv>

Known columns: every quoted word in dgva-pipeline/oracle_pipeline/submission_to_exchange/*.pl and
perl_modules/submission/*.pm (the converter reads columns by their header text). A header counts as known when its
normalised form (lower case, spaces and dashes as underscores) is one of those words.
Header row of a sheet: the first of its first 30 rows with at least 3 known headers.
Writes one row per spreadsheet, sheet and unknown column with values, with the number of non-empty cells and examples.
"""
import csv
import glob
import os
import re
import sys

import openpyxl
import xlrd

root, pipeline, paths_tsv, out = sys.argv[1:5]

known = set()
for f in glob.glob(os.path.join(pipeline, 'oracle_pipeline/submission_to_exchange/*.pl')) + \
        glob.glob(os.path.join(pipeline, 'oracle_pipeline/perl_modules/submission/*.pm')):
    known.update(w.lower() for w in re.findall(r"['\"]([A-Za-z][A-Za-z0-9_ ]{1,60})['\"]", open(f, errors='replace').read()))
    known.update(w.lower() for w in re.findall(r"\{([A-Za-z][A-Za-z0-9_]{1,60})\}", open(f, errors='replace').read()))

folder_study = {r['folder_id']: r['study_accession'] for r in csv.DictReader(open(paths_tsv), delimiter='\t') if r['folder_id']}


def norm(h):
    return re.sub(r'[\s\-]+', '_', str(h).strip()).lower().strip('_*')


def sheets(path):
    if path.lower().endswith(('.xlsx', '.xlsm')):
        wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
        for ws in wb.worksheets:
            yield ws.title, ([c for c in row] for row in ws.iter_rows(values_only=True))
    else:
        wb = xlrd.open_workbook(path)
        for ws in wb.sheets():
            yield ws.name, (ws.row_values(i) for i in range(ws.nrows))


files = [p for p in glob.glob(os.path.join(root, '**', '*.xls*'), recursive=True)
         if re.search(r'\.xlsx?$', p, re.I)]
# gzipped spreadsheets (*.xls.gz, *.xlsx.gz) are decompressed to a temporary file
files += [p for p in glob.glob(os.path.join(root, '**', '*.xls*.gz'), recursive=True)]

import gzip
import tempfile
tmpdir = tempfile.mkdtemp()

with open(out, 'w', newline='') as o:
    w = csv.writer(o, delimiter='\t')
    w.writerow(['study', 'file', 'sheet', 'header_row', 'column', 'non_empty', 'examples'])
    for p in sorted(set(files)):
        rel = os.path.relpath(p, root)
        m = re.search(r'(VLOAD_\d+|ELOAD_\d+)', rel)
        study = folder_study.get(m.group(1), '') if m else ''
        src = p
        if p.endswith('.gz'):
            src = os.path.join(tmpdir, os.path.basename(p)[:-3])
            with gzip.open(p, 'rb') as fi, open(src, 'wb') as fo:
                fo.write(fi.read())
        try:
            for sname, rows in sheets(src):
                rows = list(rows)
                hdr_i = None
                for i, row in enumerate(rows[:30]):
                    if sum(1 for c in row if c is not None and norm(c) in known) >= 3:
                        hdr_i = i
                        break
                if hdr_i is None and sum(1 for r in rows[:200] if r and r[0] is not None and norm(r[0]) in known) >= 3:
                    # vertical sheet (STUDY, SAMPLESETS): field names in column A, values to the right
                    for r in rows:
                        if not r or r[0] is None or str(r[0]).strip() == '' or norm(r[0]) in known:
                            continue
                        vals = [c for c in r[1:] if c not in (None, '')]
                        if vals:
                            ex = '; '.join(sorted({str(v)[:40] for v in vals})[:3])
                            w.writerow([study, rel, sname, 'vertical', str(r[0]).strip(), len(vals), ex])
                    continue
                if hdr_i is None:
                    if sname.startswith(('Info-', 'PLEASE READ')):
                        continue
                    if any(any(c not in (None, '') for c in r) for r in rows):
                        w.writerow([study, rel, sname, '', '(no header row recognised)', len(rows), ''])
                    continue
                hdr = rows[hdr_i]
                for j, h in enumerate(hdr):
                    if h is None or str(h).strip() == '' or norm(h) in known:
                        continue
                    vals = [r[j] for r in rows[hdr_i + 1:] if j < len(r) and r[j] not in (None, '')]
                    if vals:
                        ex = '; '.join(sorted({str(v)[:40] for v in vals})[:3])
                        w.writerow([study, rel, sname, hdr_i + 1, str(h).strip(), len(vals), ex])
        except Exception as e:
            w.writerow([study, rel, '', '', f'ERROR {type(e).__name__}: {str(e)[:80]}', '', ''])
        finally:
            if src != p and os.path.exists(src):
                os.remove(src)
print('files', len(set(files)))
