"""Match the calls of a loaded Exchange XML to the genotypes of the study's source VCF (files only, no database).

Usage: python -I match_calls_genotypes.py <loaded_exchange.xml.gz> <genotypes.vcf.gz> <out_prefix>

XML side: each VARIANT_CALL gives (linked region submitter ID, sample) and its call type and position.
VCF side: (ID, sample) -> genotype; ID is the submitter region ID.
Reports how the calls map to genotypes: carrier (het/hom/hemi), reference, missing, or no record/sample at all,
and the carriers that have no call. Writes <out_prefix>.summary.txt and <out_prefix>.examples.tsv.
"""
import collections
import gzip
import sys
import xml.etree.ElementTree as ET

xml_path, vcf_path, out = sys.argv[1:4]

calls = []
with gzip.open(xml_path, 'rb') as fh:
    for ev, el in ET.iterparse(fh, events=('end',)):
        if el.tag == 'VARIANT_CALL':
            s = el.find('SAMPLE')
            regs = [r.get('variant_region_id') for r in el.findall('LINKED_VARIANT_REGION')]
            g = el.find('PLACEMENT/GENOME')
            calls.append((el.get('variant_call_id'), regs, s.get('sample_id') if s is not None else None,
                          el.get('variant_call_type'), g.get('chr_name') if g is not None else None))
            el.clear()
        elif el.tag == 'VARIANT_REGION':
            el.clear()

needed_regions = {r for c in calls for r in c[1]}
gts = {}
vcf_samples = []
vcf_ids = set()
with gzip.open(vcf_path, 'rt') as fh:
    for line in fh:
        if line.startswith('##'):
            continue
        f = line.rstrip('\n').split('\t')
        if line.startswith('#CHROM'):
            vcf_samples = f[9:]
            continue
        vcf_ids.add(f[2])
        if f[2] not in needed_regions:
            continue
        for s, g in zip(vcf_samples, f[9:]):
            gts[(f[2], s)] = g.split(':')[0]


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


vcf_sample_set = set(vcf_samples)
res = collections.Counter()
by_type = collections.defaultdict(collections.Counter)
examples = collections.defaultdict(list)
called = set()
for cid, regs, s, ctype, chrom in calls:
    if not regs:
        k = 'call with no linked region'
    elif len(regs) > 1:
        k = 'call linked to several regions'
    elif regs[0] not in vcf_ids:
        k = 'region not in VCF'
    elif s not in vcf_sample_set:
        k = 'sample not in VCF'
    else:
        k = kind(gts[(regs[0], s)])
        called.add((regs[0], s))
    res[k] += 1
    by_type[ctype][k] += 1
    if len(examples[k]) < 5:
        examples[k].append((cid, ','.join(regs), s, ctype, chrom, gts.get((regs[0], s), '') if regs else ''))

carriers_no_call = collections.Counter()
for (r, s), g in gts.items():
    k = kind(g)
    if k in ('heterozygous', 'homozygous', 'hemizygous') and (r, s) not in called:
        carriers_no_call[k] += 1
        if len(examples['carrier with no call']) < 5:
            examples['carrier with no call'].append(('', r, s, '', '', g))

with open(out + '.summary.txt', 'w') as o:
    p = lambda *x: print(*x, file=o)
    p('xml', xml_path)
    p('vcf', vcf_path)
    p('calls:', len(calls), ' regions linked:', len(needed_regions), ' VCF records:', len(vcf_ids),
      ' VCF samples:', len(vcf_samples))
    p('calls by genotype of their (region, sample) in the VCF:')
    for k, v in res.most_common():
        p(f'  {k}: {v}')
    p('by call type:')
    for t, c in sorted(by_type.items(), key=lambda x: -sum(x[1].values())):
        p(f'  {t}: {dict(c.most_common())}')
    p('VCF carriers of the linked regions with no call:', dict(carriers_no_call))
with open(out + '.examples.tsv', 'w') as o:
    o.write('case\tcall_id\tregion\tsample\tcall_type\tchr\tvcf_gt\n')
    for k, rows in examples.items():
        for r in rows:
            o.write(k + '\t' + '\t'.join(str(x) for x in r) + '\n')
print(open(out + '.summary.txt').read())
