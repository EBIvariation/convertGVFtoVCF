"""Compare a genotype VCF (submitter region IDs in ID) with an Exchange XML's calls of the same study.

Usage: python -I compare_vcf_exchange.py <vcf.gz> <exchange.xml.gz> <out_prefix>
VCF side: (ID, sample) -> zygosity for each carrier (first ID of the record).
XML side: region variant_region_id -> its SUPPORTING_VARIANT_CALL ids; call -> sample, zygosity, and the
attributes copy_number / allele_frequency / allele_count / allele_number when present. No database access.
"""
import collections
import gzip
import sys
import xml.etree.ElementTree as ET

vcf_path, xml_path, out = sys.argv[1:4]

vcf = {}
with gzip.open(vcf_path, 'rt') as fh:
    for line in fh:
        if line.startswith('##'):
            continue
        f = line.rstrip('\n').split('\t')
        if line.startswith('#CHROM'):
            samples = f[9:]
            continue
        rid = f[2].split(';')[0]
        for s, g in zip(samples, f[9:]):
            gt = g.split(':')[0]
            al = gt.replace('|', '/').split('/')
            nz = [a for a in al if a not in ('0', '.')]
            if not nz:
                continue
            if len(al) == 1:
                z = 'hemizygous'
            elif len(nz) == 2 and nz[0] == nz[1]:
                z = 'homozygous'
            else:
                z = 'heterozygous'
            vcf[(rid, s)] = z

call_info = {}
region_calls = collections.defaultdict(list)
attr_counts = collections.Counter()
with gzip.open(xml_path, 'rb') as fh:
    cur_region = None
    for ev, el in ET.iterparse(fh, events=('start', 'end')):
        if ev == 'start':
            if el.tag == 'VARIANT_REGION':
                cur_region = el.get('variant_region_id')
            continue
        if el.tag == 'SUPPORTING_VARIANT_CALL' and cur_region:
            region_calls[cur_region].append(el.get('variant_call_id'))
        elif el.tag == 'VARIANT_CALL':
            s = el.find('SAMPLE')
            for a in ('copy_number', 'allele_frequency', 'allele_count', 'allele_number', 'zygosity',
                      'reference_copy_number'):
                if el.get(a):
                    attr_counts[a] += 1
            call_info[el.get('variant_call_id')] = (s.get('sample_id') if s is not None else None,
                                                    (el.get('zygosity') or '').lower())
            el.clear()
        elif el.tag == 'VARIANT_REGION':
            cur_region = None
            el.clear()

xml = {}
for rid, calls in region_calls.items():
    for c in calls:
        if c in call_info:
            s, z = call_info[c]
            xml[(rid, s)] = z

vk, xk = set(vcf), set(xml)
both = vk & xk
pairs = collections.Counter((vcf[k], xml[k]) for k in both)
with open(out + '.summary.txt', 'w') as o:
    p = lambda *x: print(*x, file=o)
    p('vcf', vcf_path)
    p('xml', xml_path)
    p('XML calls:', len(call_info), ' regions:', len(region_calls), ' call attributes with a value:', dict(attr_counts))
    p('VCF carriers:', len(vk), ' XML (region, sample):', len(xk))
    p('in both:', len(both), ' VCF only:', len(vk - xk), ' XML only:', len(xk - vk))
    p('zygosity pairs (VCF, XML):', dict(pairs.most_common(10)))
    p('examples VCF only:', sorted(vk - xk)[:5])
    p('examples XML only:', sorted(xk - vk)[:5])
print(open(out + '.summary.txt').read())
