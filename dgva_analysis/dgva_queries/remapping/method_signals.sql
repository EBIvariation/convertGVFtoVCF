-- Placements counted by placement method, remap score (present or not), alignment pass, placement rank and remap
-- failure code.
-- Binds: study (a study accession, or ALL). With ALL, reads every placement.
SELECT
    m.variant_placement_method,
    CASE WHEN r.variant_placement_id IS NULL THEN 'N' ELSE 'Y' END AS has_remap_score,
    a.alignment,
    pr.placement_rank,
    CASE WHEN f.variant_placement_id IS NULL THEN 'N' ELSE 'Y' END AS has_remap_failure_code,
    COUNT(*) AS placements
FROM dgva.variant_placement p
    JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
    LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_remap r ON r.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_align a ON a.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.var_plac_placement_rank pr ON pr.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.var_plac_remap_failure_code f ON f.variant_placement_id = p.variant_placement_id
WHERE :study = 'ALL' OR vs.study_accession = :study
GROUP BY m.variant_placement_method,
    CASE WHEN r.variant_placement_id IS NULL THEN 'N' ELSE 'Y' END,
    a.alignment, pr.placement_rank,
    CASE WHEN f.variant_placement_id IS NULL THEN 'N' ELSE 'Y' END
ORDER BY placements DESC
