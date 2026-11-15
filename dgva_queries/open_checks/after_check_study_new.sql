-- Ticket 8.2. VARIANT_STUDY_NEW (14.6 M rows, same shape as VARIANT_STUDY): is every row also in VARIANT_STUDY?
-- Matched by VARIANT_STUDY_ID, then compared on accession, study, type and submitter id. Per study.
-- Binds: none. Reads both tables in full: run outside busy hours.
SELECT
    n.study_accession,
    CASE WHEN v.variant_study_id IS NULL THEN 'id_not_in_variant_study'
         WHEN DECODE(n.variant_accession, v.variant_accession, 1, 0) = 1
              AND DECODE(n.study_accession, v.study_accession, 1, 0) = 1
              AND DECODE(n.region_or_call, v.region_or_call, 1, 0) = 1
              AND DECODE(n.submitter_variant_id, v.submitter_variant_id, 1, 0) = 1 THEN 'same_values'
         ELSE 'different_values' END AS status,
    CASE WHEN EXISTS (SELECT 1 FROM dgva.variant_study a WHERE a.variant_accession = n.variant_accession)
         THEN 'Y' ELSE 'N' END AS accession_in_variant_study,
    CASE WHEN w.variant_accession IS NULL THEN 'N' ELSE 'Y' END AS withdrawn,
    COUNT(*) AS variant_rows,
    MIN(n.created_date) AS first_created,
    MAX(n.created_date) AS last_created
FROM dgva.variant_study_new n
    LEFT JOIN dgva.variant_study v ON v.variant_study_id = n.variant_study_id
    LEFT JOIN dgva.variant_withdrawn w ON w.variant_accession = n.variant_accession
GROUP BY n.study_accession,
    CASE WHEN v.variant_study_id IS NULL THEN 'id_not_in_variant_study'
         WHEN DECODE(n.variant_accession, v.variant_accession, 1, 0) = 1
              AND DECODE(n.study_accession, v.study_accession, 1, 0) = 1
              AND DECODE(n.region_or_call, v.region_or_call, 1, 0) = 1
              AND DECODE(n.submitter_variant_id, v.submitter_variant_id, 1, 0) = 1 THEN 'same_values'
         ELSE 'different_values' END,
    CASE WHEN EXISTS (SELECT 1 FROM dgva.variant_study a WHERE a.variant_accession = n.variant_accession)
         THEN 'Y' ELSE 'N' END,
    CASE WHEN w.variant_accession IS NULL THEN 'N' ELSE 'Y' END
ORDER BY n.study_accession, status
