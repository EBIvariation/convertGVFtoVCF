"""Ticket 2, follow-up: look up variants with no submitted placement in a study's stored Exchange XML.

Runs locally on the outputs of the saved queries; it does not connect to the database.

1. Rebuild a STUDY_FILE from the file_chunks.sql output(s) (gzip is decompressed):

    python find_variants_in_exchange.py rebuild --chunks nstd35_file1.chunks.out --size 123456 --output nstd35_file1.xml

2. Search one or more rebuilt Exchange XML files for the variants listed by variants_without_submitted_placement.sql:

    python find_variants_in_exchange.py search --variants nstd35_variants.out --xml nstd35_file1.xml nstd35_file2.xml \
        --output nstd35_exchange_placements.tsv

The search writes one row per placement found in the XML for each listed variant (or one row with found=N), and
prints a summary: how many variants have a submitted placement in the XML, only remapped ones, or are missing.
"""
import argparse
import ast
import csv
import gzip
import sys
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict

csv.field_size_limit(sys.maxsize)

VARIANT_ELEMENTS = {
    'VARIANT_CALL': ('variant_call_accession', 'variant_call_id', 'call'),
    'VARIANT_REGION': ('variant_region_accession', 'variant_region_id', 'region'),
}
GENOME_ATTRIBUTES = ['assembly', 'chr_name', 'chr_accession', 'contig_accession', 'start', 'stop', 'outer_start',
                     'outer_stop', 'inner_start', 'inner_stop', 'remap_score', 'alignment', 'placement_rank',
                     'remap_failure_code']
NO_POSITION_METHODS = {'Submitted without placement on assembly', 'Submitted on non-INSDC assembly'}
OUTPUT_COLUMNS = ['variant_case', 'variant_accession', 'submitter_variant_id', 'region_or_call', 'db_placements',
                  'found', 'matched_on', 'xml_file', 'placement_method', 'cytoband'] + GENOME_ATTRIBUTES


def local_name(tag):
    return tag.rsplit('}', 1)[-1]


def rebuild(chunk_files, output, expected_size=None):
    chunks = {}
    for chunk_file in chunk_files:
        with open(chunk_file, newline='') as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                value = row['BYTES']
                data = ast.literal_eval(value) if value.startswith(("b'", 'b"')) else bytes.fromhex(value)
                chunks[int(row['CHUNK'])] = data
    if not chunks:
        sys.exit('No chunks read.')
    missing = sorted(set(range(1, max(chunks) + 1)) - set(chunks))
    if missing:
        sys.exit(f'Missing chunks: {missing[:10]}{" ..." if len(missing) > 10 else ""}')
    content = b''.join(chunks[n] for n in sorted(chunks))
    if expected_size is not None and len(content) != expected_size:
        sys.exit(f'Rebuilt {len(content)} bytes, expected {expected_size}. Did every chunk run complete?')
    if content[:2] == b'\x1f\x8b':
        content = gzip.decompress(content)
        print(f'Content was gzip; decompressed to {len(content)} bytes.', file=sys.stderr)
    with open(output, 'wb') as handle:
        handle.write(content)
    print(f'Wrote {output}: {len(content)} bytes from {len(chunks)} chunks.', file=sys.stderr)


def read_variants(variants_file):
    """One entry per variant accession, with the database placements summarised."""
    variants = {}
    placements = defaultdict(list)
    with open(variants_file, newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            accession = row['VARIANT_ACCESSION']
            variants[accession] = {
                'variant_case': row['VARIANT_CASE'],
                'variant_accession': accession,
                'submitter_variant_id': row['SUBMITTER_VARIANT_ID'],
                'region_or_call': row['REGION_OR_CALL'],
            }
            placements[accession].append(
                f"{row['PLACEMENT_CLASS']}:{row['ASSEMBLY_NAME'] or row['ASSEMBLY_ACCESSION']}:"
                f"{row['CHROMOSOME_NAME']}:{row['START_POS'] or row['OUTER_START_POS'] or row['INNER_START_POS']}-"
                f"{row['STOP_POS'] or row['OUTER_STOP_POS'] or row['INNER_STOP_POS']}"
            )
    for accession, variant in variants.items():
        variant['db_placements'] = ';'.join(placements[accession])
    return variants


def iter_xml_variants(xml_file):
    """Yield (region_or_call, accession, submitter id, [placement dicts]) for each variant element in the file."""
    for _, element in ET.iterparse(xml_file, events=('end',)):
        name = local_name(element.tag)
        if name not in VARIANT_ELEMENTS:
            continue
        accession_attribute, id_attribute, region_or_call = VARIANT_ELEMENTS[name]
        placements = []
        for child in element:
            if local_name(child.tag) != 'PLACEMENT':
                continue
            placement = {'placement_method': child.get('placement_method', '')}
            for sub in child:
                if local_name(sub.tag) == 'GENOME':
                    placement.update({a: sub.get(a, '') for a in GENOME_ATTRIBUTES})
                elif local_name(sub.tag) == 'CYTOGENETIC':
                    placement['cytoband'] = sub.get('cytoband', '')
            placements.append(placement)
        yield region_or_call, element.get(accession_attribute), element.get(id_attribute), placements
        element.clear()


def search(variants_file, xml_files, output):
    variants = read_variants(variants_file)
    by_submitter_id = {(v['region_or_call'], v['submitter_variant_id']): a for a, v in variants.items()}
    found = defaultdict(list)  # accession -> [(xml_file, matched_on, placement)]
    for xml_file in xml_files:
        seen = 0
        for region_or_call, accession, submitter_id, placements in iter_xml_variants(xml_file):
            seen += 1
            if accession in variants:
                key, matched_on = accession, 'accession'
            elif (region_or_call, submitter_id) in by_submitter_id:
                key, matched_on = by_submitter_id[(region_or_call, submitter_id)], 'submitter_id'
            else:
                continue
            for placement in placements or [{'placement_method': '(no PLACEMENT element)'}]:
                found[key].append((xml_file, matched_on, placement))
        print(f'{xml_file}: {seen} variant elements read.', file=sys.stderr)

    summary = Counter()
    with open(output, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter='\t', extrasaction='ignore',
                                lineterminator='\n')
        writer.writeheader()
        for accession, variant in sorted(variants.items(), key=lambda item: (item[1]['variant_case'], item[0])):
            hits = found.get(accession)
            if not hits:
                writer.writerow({**variant, 'found': 'N'})
                summary[(variant['variant_case'], 'not in the XML')] += 1
                continue
            methods = {p.get('placement_method', '') for _, _, p in hits}
            if any(m == '' or (m.startswith('Submitted') and m not in NO_POSITION_METHODS) for m in methods):
                outcome = 'submitted placement in the XML'
            elif any(m.startswith('Remapped') for m in methods):
                outcome = 'only remapped placements in the XML'
            else:
                outcome = 'only placements without position in the XML'
            summary[(variant['variant_case'], outcome)] += 1
            for xml_file, matched_on, placement in hits:
                writer.writerow({**variant, **placement, 'found': 'Y', 'matched_on': matched_on, 'xml_file': xml_file})

    print('\nVariants by case and by what the Exchange XML holds for them:', file=sys.stderr)
    for (variant_case, outcome), count in sorted(summary.items()):
        print(f'  {variant_case}\t{outcome}\t{count}', file=sys.stderr)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(dest='command', required=True)
    rebuild_parser = subparsers.add_parser('rebuild', help='Rebuild a STUDY_FILE from file_chunks.sql outputs')
    rebuild_parser.add_argument('--chunks', nargs='+', required=True, help='file_chunks.sql output file(s)')
    rebuild_parser.add_argument('--size', type=int, help='SIZE_BYTES from study_files.out, to check the result')
    rebuild_parser.add_argument('--output', required=True)
    search_parser = subparsers.add_parser('search', help='Look up the variants in rebuilt Exchange XML files')
    search_parser.add_argument('--variants', required=True, help='variants_without_submitted_placement.sql output')
    search_parser.add_argument('--xml', nargs='+', required=True, help='Rebuilt Exchange XML file(s)')
    search_parser.add_argument('--output', required=True)
    args = parser.parse_args(argv)
    if args.command == 'rebuild':
        rebuild(args.chunks, args.output, args.size)
    else:
        search(args.variants, args.xml, args.output)


if __name__ == '__main__':
    main()
