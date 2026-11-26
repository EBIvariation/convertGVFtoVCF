-- Copy numbers and allele frequencies of estd214 and estd219: rows in VARIANT_CALL_COPY_NUMBER and
-- ALLELE_FREQUENCY for each study, and non-null values in the load tables ESTD214_CALLS_CN and ESTD219_CALLS_CN.
-- Binds: none. Reads the two load tables (15.8 M rows) and the calls of the two studies.
SELECT 'archive VARIANT_CALL_COPY_NUMBER' AS source, vs.study_accession, COUNT(*) AS row_count
FROM dgva.variant_study vs JOIN dgva.variant_call_copy_number c ON c.variant_call_accession = vs.variant_accession
WHERE vs.study_accession IN ('estd214', 'estd219') GROUP BY vs.study_accession
UNION ALL
SELECT 'archive ALLELE_FREQUENCY', vs.study_accession, COUNT(*)
FROM dgva.variant_study vs JOIN dgva.allele_frequency a ON a.variant_call_accession = vs.variant_accession
WHERE vs.study_accession IN ('estd214', 'estd219') GROUP BY vs.study_accession
UNION ALL
SELECT 'load table copy_number', 'estd214', COUNT(copy_number) FROM dgva.estd214_calls_cn
UNION ALL
SELECT 'load table copy_number', 'estd219', COUNT(copy_number) FROM dgva.estd219_calls_cn
UNION ALL
SELECT 'load table af', 'estd219', COUNT(af) FROM dgva.estd219_calls_cn
