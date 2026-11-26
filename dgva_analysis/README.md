# DGVa analysis

Read-only queries against the DGVa Oracle database, a client that runs them, and scripts that read DGVa files
(GVF, Exchange XML, VCF, submission spreadsheets) without connecting to the database.

## Read-only client

`read_only_query.py` provides `ReadOnlyOracleClient`. It accepts only a single `SELECT` or `WITH` statement, rejects
DML, DDL, PL/SQL, transaction control, `SELECT ... FOR UPDATE`, database links and calls to functions other than known
Oracle built-ins, runs each statement in a read-only transaction, never commits, and refuses an autocommit
connection.

`dgva_readonly_query.py` runs one saved query and writes the result as TSV:

```
python -m dgva_analysis.dgva_readonly_query --config <config.yaml> \
    --sql dgva_analysis/dgva_queries/schema/tables.sql --bind owner=DGVA --max-rows 50 --output tables.out
```

- `--config`: YAML file with a `DGVA` section holding `host`, `port`, `user`, `password` and `service_name`.
- `--bind name=value`: a bind variable (repeatable). Queries take every value as a bind variable (`:name`).
- `--max-rows`: row limit, 1000 by default; `0` for no limit.

Tests: `python -m pytest dgva_analysis/tests`. They also check that every `.sql` file under `dgva_queries/` passes
`check_read_only`.

## Queries (`dgva_queries/`)

Each `.sql` file holds one statement. Its header comment says what it returns, its binds, and the tables it reads.

Several queries count **submitted placements**: a placement whose `VARIANT_PLACEMENT_METHOD` is null or starts with
`Submitted` (apart from `Submitted without placement on assembly` and `Submitted on non-INSDC assembly`), and that has
no `VARIANT_PLACEMENT_REMAP` row.

| Folder | Contents |
|---|---|
| `schema/` | Tables, columns, views and recorded modifications of a schema (bind `owner`) |
| `admin/` | Sessions opened by the read-only client |
| `remapping/` | Placement methods and the other signals of a remapping; placement classes per study and assembly; variants with no submitted placement; DGVa's remapping tables |
| `remapping/exchange_xml/` | Files stored in `STUDY_FILE`, read in 2000-byte chunks, and the variants to look up in them |
| `whole_variants/` | Submitted placements flagged with the conditions under which the DGVa GVF exporter writes them; region-to-call links; per-study accession lists |
| `open_checks/` | Confidence-interval tables compared with each other and with inner and outer coordinates; one-off load and withdrawal tables compared with the archive tables; per-study counts of fields absent from GVF files |
| `gap_table/` | Per-study counts of values in variant-level, placement-level and metadata columns; `checks/` lists the values of individual columns |
| `load_loss/` | Load records (`STAGE_LOG`, `STAGE_TO_ARCHIVE_LOG`, `STAGE_TO_ARCHIVE_ROWCOUNT`), current row counts of the call tables, and variants present only in `VARIANT_STUDY_NEW` |
| `report/` | Short queries: zygosity per study, chromosome X zygosity, row counts of empty tables, copy numbers and allele frequencies, duplicated allele rows, values from submission spreadsheets |
| `storage/` | VLOAD ID and display name of every study |

## Scripts that read files

None of these connects to the database. Each one's docstring gives its arguments and output.

| Script | Reads | Writes |
|---|---|---|
| `remapping/exchange_xml/find_variants_in_exchange.py` | Output of `file_chunks.sql` and `variants_without_submitted_placement.sql` | The rebuilt file; the placements of the listed variants found in it |
| `whole_variants/gvf_counts.py` | GVF files | Lines and distinct accessions per file, for calls and regions |
| `open_checks/gvf_bug_checks.py` | GVF files | Header keys and column 9 attributes per file, with derived checks |
| `load_loss/scan_exchange.py` | Exchange XML files | Element and attribute counts per study |
| `load_loss/marker_search.py` | Exchange XML files | Counts of genotype, inheritance, event, mate and group ID markers per file |
| `load_loss/compare_genotypes.py` | A genotype VCF and a GVF | Carrier genotypes compared with GVF calls by region and sample |
| `load_loss/compare_vcf_exchange.py` | A genotype VCF and an Exchange XML | Carrier genotypes compared with the XML's calls by region and sample |
| `load_loss/match_calls_genotypes.py` | An Exchange XML and a genotype VCF | Each call classified by its sample's genotype in its region |
| `load_loss/match_calls_genotypes_linked.py` | An Exchange XML and genotype VCFs | The same, through region-to-region links |
| `load_loss/spreadsheet_columns.py` | Submission spreadsheets and the `dgva-pipeline` converter code | Spreadsheet columns the converter never names, with their non-empty cells |

`spreadsheet_columns.py` needs `openpyxl` and `xlrd`.

Query and script outputs (`*.out`, `*.tsv`, `*.log`, `*.summary.txt`) and files rebuilt from `STUDY_FILE` are not
tracked (`dgva_queries/.gitignore`).
