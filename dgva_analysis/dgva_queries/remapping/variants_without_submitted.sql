-- Per study: non-withdrawn variants counted by the classes of placements they have (see study_assembly_class.sql;
-- a "conflict" placement counts as neither submitted nor remapped).
-- null_and_named_submitted: variants with both a null-method placement and a "Submitted ..." placement.
-- Binds: study (a study accession, or ALL). With ALL, reads every placement.
SELECT
    study_accession,
    region_or_call,
    SUM(CASE WHEN n_submitted > 0 THEN 1 ELSE 0 END) AS with_submitted,
    SUM(CASE WHEN n_submitted = 0 AND n_remapped > 0 THEN 1 ELSE 0 END) AS remapped_only,
    SUM(CASE WHEN n_submitted = 0 AND n_remapped = 0 AND n_placements > 0 THEN 1 ELSE 0 END) AS other_placements_only,
    SUM(CASE WHEN n_placements = 0 THEN 1 ELSE 0 END) AS no_placement,
    SUM(CASE WHEN n_submitted_assemblies > 1 THEN 1 ELSE 0 END) AS submitted_on_several_assemblies,
    SUM(CASE WHEN n_null_method > 0 AND n_named_submitted > 0 THEN 1 ELSE 0 END) AS null_and_named_submitted,
    COUNT(*) AS variants
FROM (
    SELECT
        vs.variant_study_id,
        vs.study_accession,
        vs.region_or_call,
        COUNT(p.variant_placement_id) AS n_placements,
        SUM(CASE WHEN p.variant_placement_id IS NOT NULL AND r.variant_placement_id IS NULL
                  AND (m.variant_placement_method IS NULL
                       OR (m.variant_placement_method LIKE 'Submitted%'
                           AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                                  'Submitted on non-INSDC assembly')))
                 THEN 1 ELSE 0 END) AS n_submitted,
        SUM(CASE WHEN m.variant_placement_method LIKE 'Remapped%' THEN 1 ELSE 0 END) AS n_remapped,
        COUNT(DISTINCT CASE WHEN p.variant_placement_id IS NOT NULL AND r.variant_placement_id IS NULL
                  AND (m.variant_placement_method IS NULL
                       OR (m.variant_placement_method LIKE 'Submitted%'
                           AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                                  'Submitted on non-INSDC assembly')))
                 THEN p.assembly_accession END) AS n_submitted_assemblies,
        SUM(CASE WHEN p.variant_placement_id IS NOT NULL AND m.variant_placement_method IS NULL
                  AND r.variant_placement_id IS NULL THEN 1 ELSE 0 END) AS n_null_method,
        SUM(CASE WHEN m.variant_placement_method LIKE 'Submitted%' THEN 1 ELSE 0 END) AS n_named_submitted
    FROM dgva.variant_study vs
        LEFT JOIN dgva.variant_placement p ON p.variant_study_id = vs.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.variant_placement_remap r ON r.variant_placement_id = p.variant_placement_id
    WHERE (:study = 'ALL' OR vs.study_accession = :study)
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
    GROUP BY vs.variant_study_id, vs.study_accession, vs.region_or_call
)
GROUP BY study_accession, region_or_call
ORDER BY study_accession, region_or_call
