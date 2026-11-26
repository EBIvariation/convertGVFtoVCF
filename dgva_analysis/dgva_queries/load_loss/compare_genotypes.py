"""Compare a genotype VCF whose ID column holds region accessions with the calls of a DGVa GVF of the same study.

Usage: python -I compare_genotypes.py <vcf.gz> <gvf.gz> <out_prefix>
VCF side: one (esv, sample, zygosity) per non-reference allele carried by a sample.
GVF side: one (parent esv, sample_name, Zygosity) per distinct call accession (Name).
Writes <out_prefix>.summary.txt and example lines; no database access.
"""
import collections
import gzip
import sys

vcf_path, gvf_path, out = sys.argv[1:4]
FIRST_ID = len(sys.argv) > 4 and sys.argv[4] == "first_id"

vcf = {}
phased = unphased = 0
fmt_keys = collections.Counter()
per_sample_cn = 0
with gzip.open(vcf_path, 'rt') as fh:
    samples = None
    for line in fh:
        if line.startswith('##'):
            continue
        f = line.rstrip('\n').split('\t')
        if line.startswith('#CHROM'):
            samples = f[9:]
            continue
        ids = f[2].split(';')
        fmt = f[8].split(':')
        fmt_keys[f[8]] += 1
        cn_i = fmt.index('CN') if 'CN' in fmt else None
        for s, g in zip(samples, f[9:]):
            parts = g.split(':')
            gt = parts[0]
            sep = '|' if '|' in gt else '/'
            al = [a for a in gt.split(sep)]
            if cn_i is not None and len(parts) > cn_i and parts[cn_i] not in ('.', ''):
                per_sample_cn += 1
            nz = [a for a in al if a not in ('0', '.')]
            if not nz:
                continue
            if sep == '|':
                phased += 1
            else:
                unphased += 1
            for a in set(nz):
                k = int(a)
                acc = ids[0] if FIRST_ID else (ids[k - 1] if len(ids) > 1 else ids[0])
                zyg = 'homozygous' if al.count(a) == 2 or (len(set(nz)) == 1 and len(nz) == 2) else ('hemizygous' if len(al) == 1 else 'heterozygous')
                vcf[(acc, s)] = zyg

gvf = {}
with gzip.open(gvf_path, 'rt') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        a = dict(kv.split('=', 1) for kv in f[8].rstrip(';').split(';') if '=' in kv)
        if 'parent' not in a or 'sample_name' not in a:
            continue
        gvf[a['Name']] = (a['parent'], a['sample_name'], a.get('Zygosity', ''))

gvf_keys = collections.Counter((p, s) for p, s, z in gvf.values())
gvf_by_key = {(p, s): z for p, s, z in gvf.values()}
vk, gk = set(vcf), set(gvf_by_key)
both = vk & gk
zyg_agree = sum(1 for k in both if vcf[k] == gvf_by_key[k])
zyg_pairs = collections.Counter((vcf[k], gvf_by_key[k]) for k in both)
vcf_accs = {a for a, s in vk}
gvf_accs = {a for a, s in gk}

with open(out + '.summary.txt', 'w') as o:
    p = lambda *x: print(*x, file=o)
    p('vcf', vcf_path)
    p('gvf', gvf_path)
    p('VCF FORMAT columns:', dict(fmt_keys))
    p('VCF carrier genotypes (sample x accession):', len(vk), ' phased:', phased, ' unphased:', unphased)
    p('VCF per-sample CN values (any sample, incl. ref):', per_sample_cn)
    p('GVF distinct calls:', len(gvf), ' distinct (parent, sample):', len(gk),
      ' (parent, sample) with >1 call:', sum(1 for v in gvf_keys.values() if v > 1))
    p('accessions: VCF', len(vcf_accs), ' GVF', len(gvf_accs), ' both', len(vcf_accs & gvf_accs))
    p('(accession, sample) in both:', len(both), ' VCF only:', len(vk - gk), ' GVF only:', len(gk - vk))
    p('zygosity agrees on', zyg_agree, 'of', len(both))
    p('zygosity pairs (VCF, GVF):', dict(zyg_pairs.most_common(10)))
    p('examples VCF only:', sorted(vk - gk)[:5])
    p('examples GVF only:', sorted(gk - vk)[:5])
print(open(out + '.summary.txt').read())
