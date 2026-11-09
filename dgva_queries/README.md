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
