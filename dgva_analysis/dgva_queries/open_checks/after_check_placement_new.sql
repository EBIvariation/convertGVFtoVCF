-- VARIANT_PLACEMENT_NEW rows classified against VARIANT_PLACEMENT: identical, different values, or missing. For
-- missing rows: whether the variant is still in VARIANT_STUDY, whether it is withdrawn, and whether the placement's
-- coordinates are in VARIANT_PLACEMENT_START.
-- Binds: none. Reads VARIANT_PLACEMENT_NEW and VARIANT_PLACEMENT in full.
SELECT
    CASE WHEN p.variant_placement_id IS NULL THEN 'not_in_variant_placement'
         WHEN DECODE(n.chromosome_name, p.chromosome_name, 1, 0) = 1
              AND DECODE(n.assembly_accession, p.assembly_accession, 1, 0) = 1
              AND DECODE(n.variant_study_id, p.variant_study_id, 1, 0) = 1
              AND DECODE(n.chr_accession, p.chr_accession, 1, 0) = 1 THEN 'same_values'
         ELSE 'different_values' END AS status,
    CASE WHEN vs.variant_study_id IS NULL THEN 'N' ELSE 'Y' END AS variant_in_variant_study,
    CASE WHEN w.variant_accession IS NULL THEN 'N' ELSE 'Y' END AS variant_withdrawn,
    CASE WHEN st.variant_placement_id IS NULL THEN 'N' ELSE 'Y' END AS has_start_row,
    COUNT(*) AS placement_rows,
    COUNT(DISTINCT n.variant_study_id) AS variants,
    MIN(vs.study_accession) AS first_study,
    MAX(vs.study_accession) AS last_study
FROM dgva.variant_placement_new n
    LEFT JOIN dgva.variant_placement p ON p.variant_placement_id = n.variant_placement_id
    LEFT JOIN dgva.variant_study vs ON vs.variant_study_id = n.variant_study_id
    LEFT JOIN dgva.variant_withdrawn w ON w.variant_accession = vs.variant_accession
    LEFT JOIN dgva.variant_placement_start st ON st.variant_placement_id = n.variant_placement_id
GROUP BY
    CASE WHEN p.variant_placement_id IS NULL THEN 'not_in_variant_placement'
         WHEN DECODE(n.chromosome_name, p.chromosome_name, 1, 0) = 1
              AND DECODE(n.assembly_accession, p.assembly_accession, 1, 0) = 1
              AND DECODE(n.variant_study_id, p.variant_study_id, 1, 0) = 1
              AND DECODE(n.chr_accession, p.chr_accession, 1, 0) = 1 THEN 'same_values'
         ELSE 'different_values' END,
    CASE WHEN vs.variant_study_id IS NULL THEN 'N' ELSE 'Y' END,
    CASE WHEN w.variant_accession IS NULL THEN 'N' ELSE 'Y' END,
    CASE WHEN st.variant_placement_id IS NULL THEN 'N' ELSE 'Y' END
ORDER BY placement_rows DESC
