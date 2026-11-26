-- Per study and zygosity, the calls of estd214 and estd219 with a placement on chromosome X that has no
-- VARIANT_PLACEMENT_REMAP row.
-- Binds: none.
SELECT vs.study_accession, z.zygosity, COUNT(DISTINCT vs.variant_accession) AS calls
FROM dgva.variant_study vs
    JOIN dgva.variant_placement p ON p.variant_study_id = vs.variant_study_id
    JOIN dgva.variant_call_zygosity z ON z.variant_call_accession = vs.variant_accession
WHERE vs.study_accession IN ('estd214', 'estd219')
    AND vs.region_or_call = 'call'
    AND p.chromosome_name = 'X'
    AND NOT EXISTS (SELECT 1 FROM dgva.variant_placement_remap r WHERE r.variant_placement_id = p.variant_placement_id)
GROUP BY vs.study_accession, z.zygosity
ORDER BY 1, 2
