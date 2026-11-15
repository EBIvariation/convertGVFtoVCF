-- Ticket 7. Check of the instant 0-row result of gap_placement_table_counts.sql for VARIANT_PLACEMENT.CHR_ACCESSION:
-- how many placements have a value at all, split by submitted / remapped (ticket 2 rule, simplified: remap row exists).
-- No binds. Reads VARIANT_PLACEMENT once.
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
