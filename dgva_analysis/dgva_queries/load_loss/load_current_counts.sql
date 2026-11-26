-- Rows per study in VARIANT_CALL, VARIANT_CALL_SAMPLE, VARIANT_CALL_SAMPLESET, VARIANT_CALL_ZYGOSITY,
-- VARIANT_CALL_COPY_NUMBER, ALLELE_FREQUENCY, ALLELE_COUNT and ALLELE_NUMBER, counted over all calls of the study,
-- withdrawn ones included. Comparable with load_rowcounts.sql.
-- Binds: study (one study accession, or ALL). Reads VARIANT_CALL_SAMPLE (about 30 M rows) and the other call tables.
WITH c AS (
    SELECT vs.study_accession, vs.variant_accession
    FROM dgva.variant_study vs
    WHERE vs.region_or_call = 'call'
        AND (:study = 'ALL' OR vs.study_accession = :study)
)
SELECT c.study_accession, 'variant_call' AS table_name, COUNT(*) AS current_rows
FROM c GROUP BY c.study_accession
UNION ALL
SELECT c.study_accession, 'variant_call_sample', COUNT(*)
FROM c JOIN dgva.variant_call_sample t ON t.variant_call_accession = c.variant_accession GROUP BY c.study_accession
UNION ALL
SELECT c.study_accession, 'variant_call_sampleset', COUNT(*)
FROM c JOIN dgva.variant_call_sampleset t ON t.variant_call_accession = c.variant_accession GROUP BY c.study_accession
UNION ALL
SELECT c.study_accession, 'variant_call_zygosity', COUNT(*)
FROM c JOIN dgva.variant_call_zygosity t ON t.variant_call_accession = c.variant_accession GROUP BY c.study_accession
UNION ALL
SELECT c.study_accession, 'variant_call_copy_number', COUNT(*)
FROM c JOIN dgva.variant_call_copy_number t ON t.variant_call_accession = c.variant_accession GROUP BY c.study_accession
UNION ALL
SELECT c.study_accession, 'allele_frequency', COUNT(*)
FROM c JOIN dgva.allele_frequency t ON t.variant_call_accession = c.variant_accession GROUP BY c.study_accession
UNION ALL
SELECT c.study_accession, 'allele_count', COUNT(*)
FROM c JOIN dgva.allele_count t ON t.variant_call_accession = c.variant_accession GROUP BY c.study_accession
UNION ALL
SELECT c.study_accession, 'allele_number', COUNT(*)
FROM c JOIN dgva.allele_number t ON t.variant_call_accession = c.variant_accession GROUP BY c.study_accession
ORDER BY 1, 2
