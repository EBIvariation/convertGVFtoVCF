"""Match calls to genotypes through the regions they are linked to, including the regions above them through
region-to-region links. Reads files only.

Usage: python -I match_calls_genotypes_linked.py <exchange.xml.gz> <out_prefix> <genotypes.vcf.gz> [<genotypes.vcf.gz> ...]

XML: VARIANT_CALL (accession, submitter id, sample); VARIANT_REGION (submitter id, LINKED_VARIANT_CALL accessions,
LINKED_VARIANT_REGION). Region links name the child region (by accession when the file has accessions); each call is
matched through its region and all the regions above it.
VCF: (ID, sample) -> genotype, for every VCF given. IDs joined with '&' are split; IDs are compared case-insensitively.
For each call: the genotypes of its sample in every genotyped region it belongs to.
"""
import collections
import gzip
import sys
import xml.etree.ElementTree as ET

xml_path, out = sys.argv[1:3]
vcf_paths = sys.argv[3:]

calls = {}                                  # accession -> (submitter id, sample)
region_calls = collections.defaultdict(set)  # region id -> call accessions
region_parents = collections.defaultdict(set)  # region key -> parent region keys
region_name = {}                               # region key (accession, else submitter id) -> submitter id
with gzip.open(xml_path, 'rb') as fh:
    for ev, el in ET.iterparse(fh, events=('end',)):
        if el.tag == 'VARIANT_CALL':
            s = el.find('SAMPLE')
            calls[el.get('variant_call_accession') or el.get('variant_call_id')] = (
                el.get('variant_call_id'), s.get('sample_id') if s is not None else None)
            el.clear()
        elif el.tag == 'VARIANT_REGION':
            rid = el.get('variant_region_accession') or el.get('variant_region_id')
            region_name[rid] = el.get('variant_region_id')
            for lc in el.findall('LINKED_VARIANT_CALL'):
                region_calls[rid].add(lc.get('variant_call_id'))
            for lr in el.findall('LINKED_VARIANT_REGION'):
                region_parents[lr.get('variant_region_id')].add(rid)
            el.clear()

gts = {}
conflicts = 0
vcf_samples_all = set()
vcf_ids = set()
for path in vcf_paths:
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('##'):
                continue
            f = line.rstrip('\n').split('\t')
            if line.startswith('#CHROM'):
                samples = f[9:]
                vcf_samples_all.update(samples)
                continue
            for vid in f[2].split('&'):
                vid = vid.upper()
                vcf_ids.add(vid)
                for s, g in zip(samples, f[9:]):
                    gt = g.split(':')[0]
                    old = gts.get((vid, s))
                    if old is not None and old.replace('|', '/') != gt.replace('|', '/'):
                        conflicts += 1
                    gts.setdefault((vid, s), gt)


def kind(gt):
    al = gt.replace('|', '/').split('/')
    nz = [a for a in al if a not in ('0', '.')]
    if all(a == '.' for a in al):
        return 'missing'
    if not nz:
        return 'reference'
    if len(al) == 1:
        return 'hemizygous'
    return 'homozygous' if len(nz) == 2 and nz[0] == nz[1] else 'heterozygous'


# regions of each call, including the regions above them
call_regions = collections.defaultdict(set)
for rid, accs in region_calls.items():
    stack, seen = [rid], set()
    while stack:
        r = stack.pop()
        if r in seen:
            continue
        seen.add(r)
        stack.extend(region_parents.get(r, ()))
    for a in accs:
        call_regions[a].update(seen)

res = collections.Counter()
by_set = collections.defaultdict(collections.Counter)
examples = collections.defaultdict(list)
for acc, (cid, s) in calls.items():
    regs = call_regions.get(acc, set())
    regs = {region_name.get(r, r) for r in regs}
    gen = [r for r in regs if r.upper() in vcf_ids]
    cset = cid.rsplit('_', 1)[0] if cid else ''
    cset = ''.join('N' if ch.isdigit() else ch for ch in cset)
    while 'NN' in cset:
        cset = cset.replace('NN', 'N')
    if not regs:
        k = 'call in no region'
    elif not gen:
        k = 'no genotyped region'
    elif s not in vcf_samples_all:
        k = 'sample not in VCFs'
    else:
        kinds = {kind(gts[(r.upper(), s)]) if (r.upper(), s) in gts else 'sample not genotyped for region'
                 for r in gen}
        carrier = kinds & {'heterozygous', 'homozygous', 'hemizygous'}
        if len(carrier) == 1:
            k = carrier.pop()
        elif len(carrier) > 1:
            k = 'carrier, several zygosities'
        else:
            k = sorted(kinds)[0]
    res[k] += 1
    by_set[cset][k] += 1
    if len(examples[k]) < 5:
        examples[k].append((acc, cid, s, ','.join(sorted(gen)) or ','.join(sorted(regs))[:80]))

with open(out + '.summary.txt', 'w') as o:
    p = lambda *x: print(*x, file=o)
    p('xml', xml_path)
    for v in vcf_paths:
        p('vcf', v)
    p('calls:', len(calls), ' regions:', len(set(region_calls) | set(region_parents)),
      ' VCF IDs:', len(vcf_ids), ' VCF samples:', len(vcf_samples_all), ' conflicting duplicate genotypes:', conflicts)
    p('calls by their genotype:')
    for k, v in res.most_common():
        p(f'  {k}: {v}')
    p('by call set (submitter ID pattern):')
    for t, c in sorted(by_set.items(), key=lambda x: -sum(x[1].values())):
        p(f'  {t}: {dict(c.most_common())}')
with open(out + '.examples.tsv', 'w') as o:
    o.write('case\tcall_accession\tcall_id\tsample\tregions\n')
    for k, rows in examples.items():
        for r in rows:
            o.write(k + '\t' + '\t'.join(str(x) for x in r) + '\n')
print(open(out + '.summary.txt').read())
