-- Per study, for nstd126, nstd134, nstd137 and nstd144: ALLELE_FREQUENCY rows, calls with a row, and distinct
-- (call, value) pairs.
-- Binds: none.
SELECT vs.study_accession,
       COUNT(*) AS rows_in_table,
       COUNT(DISTINCT a.variant_call_accession) AS calls,
       COUNT(DISTINCT a.variant_call_accession || '|' || a.allele_frequency) AS distinct_call_values
FROM dgva.variant_study vs
    JOIN dgva.allele_frequency a ON a.variant_call_accession = vs.variant_accession
WHERE vs.study_accession IN ('nstd126', 'nstd134', 'nstd137', 'nstd144')
GROUP BY vs.study_accession
ORDER BY 1
