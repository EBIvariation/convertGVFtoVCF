-- DGVA_REMAPPED_PLACEMENTS rows by placement method, and whether the placement they point to is in
-- VARIANT_PLACEMENT.
-- Binds: none. Reads about 40 M rows.
SELECT
    d.variant_placement_method,
    CASE WHEN p.variant_placement_id IS NULL THEN 'N' ELSE 'Y' END AS placement_still_exists,
    COUNT(*) AS placements
FROM dgva.dgva_remapped_placements d
    LEFT JOIN dgva.variant_placement p ON p.variant_placement_id = d.variant_placement_id
GROUP BY d.variant_placement_method, CASE WHEN p.variant_placement_id IS NULL THEN 'N' ELSE 'Y' END
ORDER BY placements DESC
