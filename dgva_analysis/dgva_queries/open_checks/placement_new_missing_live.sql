-- VARIANT_PLACEMENT_NEW rows missing from VARIANT_PLACEMENT whose variant is in VARIANT_STUDY, per study, variant
-- kind and assembly, with whether the variant has a placement on the same assembly in VARIANT_PLACEMENT.
-- Binds: none.
SELECT
    vs.study_accession,
    vs.region_or_call,
    n.assembly_accession,
    ga.assembly_name,
    COUNT(*) AS placement_rows,
    SUM(CASE WHEN EXISTS (SELECT 1 FROM dgva.variant_placement p WHERE p.variant_study_id = n.variant_study_id
                              AND p.assembly_accession = n.assembly_accession) THEN 1 ELSE 0 END)
        AS variant_has_placement_on_same_assembly
FROM dgva.variant_placement_new n
    JOIN dgva.variant_study vs ON vs.variant_study_id = n.variant_study_id
    LEFT JOIN dgva.genome_assembly ga ON ga.assembly_accession = n.assembly_accession
WHERE NOT EXISTS (SELECT 1 FROM dgva.variant_placement p WHERE p.variant_placement_id = n.variant_placement_id)
GROUP BY vs.study_accession, vs.region_or_call, n.assembly_accession, ga.assembly_name
ORDER BY vs.study_accession
