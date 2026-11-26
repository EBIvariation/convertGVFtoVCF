-- For each of the tables ESTD214_REGIONS_TO_WITHDRAW, ESTD222_WITHDRAW and ESTD225_WITHDRAW_REGIONS: rows, distinct
-- accessions, how many are in VARIANT_WITHDRAWN and in VARIANT_STUDY, how many are in VARIANT_STUDY but not withdrawn,
-- and the withdrawal reasons.
-- Binds: none. Small tables only.
WITH lists AS (
    SELECT 'ESTD214_REGIONS_TO_WITHDRAW' AS list_table, variant_region_accession AS variant_accession
    FROM dgva.estd214_regions_to_withdraw
    UNION ALL
    SELECT 'ESTD222_WITHDRAW', variant_accession FROM dgva.estd222_withdraw
    UNION ALL
    SELECT 'ESTD225_WITHDRAW_REGIONS', variant_accession FROM dgva.estd225_withdraw_regions
)
SELECT
    l.list_table,
    COUNT(*) AS list_rows,
    COUNT(DISTINCT l.variant_accession) AS distinct_accessions,
    SUM(CASE WHEN w.variant_accession IS NOT NULL THEN 1 ELSE 0 END) AS in_variant_withdrawn,
    SUM(CASE WHEN vs.variant_accession IS NOT NULL THEN 1 ELSE 0 END) AS in_variant_study,
    SUM(CASE WHEN w.variant_accession IS NULL AND vs.variant_accession IS NOT NULL THEN 1 ELSE 0 END)
        AS still_live_not_withdrawn,
    LISTAGG(DISTINCT w.withdrawal_reason, ' | ') WITHIN GROUP (ORDER BY w.withdrawal_reason) AS reasons
FROM lists l
    LEFT JOIN dgva.variant_withdrawn w ON w.variant_accession = l.variant_accession
    LEFT JOIN dgva.variant_study vs ON vs.variant_accession = l.variant_accession
GROUP BY l.list_table
ORDER BY l.list_table
