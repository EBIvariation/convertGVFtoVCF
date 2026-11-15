# DGVa queries

Saved, reviewed queries for the DGVa gap analysis (EVA-2269). Each file holds one `SELECT` or `WITH` statement and is
run through the read-only helper (`convert_gvf_to_vcf/metadata_retrievers/read_only_query.py`), for example:

```
dgva_readonly_query.py --config <config outside the repo> --sql dgva_queries/schema/tables.sql --bind owner=DGVA
```

Rules (from the sprint tickets):

- Values go in bind variables (`:name`), never in the SQL text.
- A second developer reviews a query before it is run against large studies.
- Keep the default row limit while developing a query. Run long queries outside busy hours.
- The same files are meant to be rerun on the cold-storage copy later.

## `remapping/exchange_xml/`: variants with no submitted placement, looked up in the study's Exchange XML

For a study such as `nstd35`, where 976 variants have only remapped placements (ticket 2). Run from
`dgva_queries/remapping/exchange_xml/`, with `--max-rows 0`:

1. `study_files.sql --bind study=nstd35`: the stored files, their size and how they are stored (`FIRST_BYTES`).
   The Exchange XMLs have type `dbvar_exchange_xml`.
2. `variants_without_submitted_placement.sql --bind study=nstd35`: every placement of the variants with no submitted
   placement (976 remapped-only and 590 with no position for `nstd35`).
3. For each `dbvar_exchange_xml` file, `file_chunks.sql --bind study=nstd35 --bind file_name=<FILE_NAME>
   --bind first_chunk=1 --bind last_chunk=<CHUNKS>`. Split a very large file into several runs.
4. Locally, with no database connection:
   `python find_variants_in_exchange.py rebuild --chunks <chunks.out> --size <SIZE_BYTES> --output <file>.xml`, then
   `python find_variants_in_exchange.py search --variants <step 2 .out> --xml <file1>.xml <file2>.xml --output
   nstd35_exchange_placements.tsv`. The summary says how many variants have a submitted placement in the XML.

## `whole_variants/`: what the GVF exporter drops (sprint 2, ticket 6)

- `exporter_drops.sql` (`--bind study=ALL`, about 30 minutes): the exporter's W1–W3 conditions applied to submitted
  placements, per study, assembly and variant type.
- `linked_variant_kinds.sql`: `LINKED_VARIANT` rows by kind (region to call, region to region) and study.
- `variant_accessions.sql`, `withdrawn_accessions.sql` (one study each): per-variant lists, compared with the FTP files.
- `non_human_sample_data.sql`: non-human studies with sample, zygosity and copy-number data.
- `gvf_counts.py` runs locally on downloaded GVFs, with no database connection. It counts distinct accessions per file
  and can write them out (`--names-dir`) for the comparison.

## `open_checks/`: confidence intervals, "after check" tables, bugs on real files (sprint 2, ticket 8)

- `ci_compare.sql` compares the two CI sets and checks them against inner and outer coordinates. `ci_examples.sql`
  and `ci_signed_check.sql` work on one study.
- The `after_check_*.sql` queries check the one-off tables against the archive tables. `placement_new_missing_live.sql`
  and `estd214_regions_to_withdraw_live.sql` explain the remainders.
- `sample_field_presence.sql` (one study) counts the database values for the fields the B bugs drop.
- `gvf_bug_checks.py` runs locally on downloaded GVFs and tallies header keys and column 9 attributes per file.

## `gap_table/`: counts for the field-level gap table (sprint 2, ticket 7)

All queries count over all studies, on submitted placements of non-withdrawn variants (ticket 2 rule).
Run them outside busy hours.

- `gap_metadata_counts.sql` (no binds, seconds): metadata candidate columns, rows with a value per study.
- `gap_variant_counts.sql` (`--bind study=ALL` or one study, about 14 minutes): variant-level columns per study and
  variant type.
- `gap_placement_table_counts.sql`: placement-level columns, one table per run. Use `--bind table_name=<TABLE>
  --bind study=ALL`; most runs take seconds to 7 minutes, and all 19 take about 45 minutes.
- `gap_placement_counts.sql`: superseded. Its single scan didn't finish in 10 hours; do not run it.
- `chr_accession_check.sql` (no binds): confirms that `VARIANT_PLACEMENT.CHR_ACCESSION` is empty everywhere.
- Outputs: `*_ALL.out`, `gap_placement_table_counts_<TABLE>.out`, and `*_nstd102.out` from the one-study test.
- `checks/` settles the checks the first gap table left open. Each query takes seconds to 2 minutes and has no binds:
  `submitted_vs_exchange_type.sql`, `alt_status_values.sql`, `linked_object_type_values.sql`,
  `assertion_method_id.sql`, `ci_in_sign.sql` and `region_parents.sql`.
