-- Placements with a non-null VARIANT_PLACEMENT.CHR_ACCESSION, split by whether the placement has a
-- VARIANT_PLACEMENT_REMAP row, with the number of variants.
-- Binds: none. Reads VARIANT_PLACEMENT once.
SELECT CASE WHEN EXISTS (SELECT 1 FROM dgva.variant_placement_remap rm
                         WHERE rm.variant_placement_id = p.variant_placement_id)
            THEN 'has remap row' ELSE 'no remap row' END AS placement_kind,
    COUNT(*) AS placements_with_chr_accession,
    COUNT(DISTINCT p.variant_study_id) AS variants
FROM dgva.variant_placement p
WHERE p.chr_accession IS NOT NULL
GROUP BY CASE WHEN EXISTS (SELECT 1 FROM dgva.variant_placement_remap rm
                           WHERE rm.variant_placement_id = p.variant_placement_id)
              THEN 'has remap row' ELSE 'no remap row' END
