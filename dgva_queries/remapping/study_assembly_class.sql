-- Ticket 2. For each study and assembly: placements in each class of the proposed rule, and how many of them
-- belong to withdrawn variants (the exporters leave those out).
--   submitted:    placement method is null or starts with "Submitted", and there is no remap score
--   remapped:     placement method starts with "Remapped"
--   no_position:  "Remap failed", "Submitted without placement on assembly", "Submitted on non-INSDC assembly"
--   conflict:     submitted by method, but carries a remap score
-- Binds: study (a study accession, or ALL). With ALL it reads all placements: run outside busy hours.
SELECT
    vs.study_accession,
    ga.assembly_name,
    p.assembly_accession,
    SUM(CASE WHEN cls = 'submitted' THEN 1 ELSE 0 END) AS submitted,
    SUM(CASE WHEN cls = 'remapped' THEN 1 ELSE 0 END) AS remapped,
    SUM(CASE WHEN cls = 'no_position' THEN 1 ELSE 0 END) AS no_position,
    SUM(CASE WHEN cls = 'conflict' THEN 1 ELSE 0 END) AS conflict,
    SUM(CASE WHEN cls = 'other' THEN 1 ELSE 0 END) AS other_method,
    SUM(CASE WHEN method IS NULL THEN 1 ELSE 0 END) AS method_null,
    SUM(CASE WHEN w.variant_accession IS NOT NULL THEN 1 ELSE 0 END) AS withdrawn,
    COUNT(*) AS placements
FROM (
    SELECT
        p.variant_placement_id,
        p.variant_study_id,
        p.assembly_accession,
        m.variant_placement_method AS method,
        CASE
            WHEN m.variant_placement_method IN ('Remap failed', 'Submitted without placement on assembly',
                                                'Submitted on non-INSDC assembly') THEN 'no_position'
            WHEN m.variant_placement_method LIKE 'Remapped%' THEN 'remapped'
            WHEN (m.variant_placement_method IS NULL OR m.variant_placement_method LIKE 'Submitted%')
                 AND r.variant_placement_id IS NOT NULL THEN 'conflict'
            WHEN m.variant_placement_method IS NULL OR m.variant_placement_method LIKE 'Submitted%' THEN 'submitted'
            ELSE 'other'
        END AS cls
    FROM dgva.variant_placement p
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.variant_placement_remap r ON r.variant_placement_id = p.variant_placement_id
) p
    JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
    LEFT JOIN dgva.genome_assembly ga ON ga.assembly_accession = p.assembly_accession
    LEFT JOIN dgva.variant_withdrawn w ON w.variant_accession = vs.variant_accession
WHERE :study = 'ALL' OR vs.study_accession = :study
GROUP BY vs.study_accession, ga.assembly_name, p.assembly_accession
ORDER BY vs.study_accession, ga.assembly_name
