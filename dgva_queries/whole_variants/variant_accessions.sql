-- Ticket 6, FTP check. One row per placement of a study's variants: accession, type, assembly, the ticket 2 class
-- of the placement, whether the variant is withdrawn, and the exporter conditions W1 (call linked from a region of
-- the same study) and W3 (start, chromosome or contig name). Compared with the accessions in the FTP GVFs.
-- Binds: study (one study accession; the output has one row per placement).
SELECT
    vs.variant_accession,
    vs.region_or_call,
    ga.assembly_name,
    p.assembly_accession,
    CASE
        WHEN p.variant_placement_id IS NULL THEN 'none'
        WHEN m.variant_placement_method IN ('Remap failed', 'Submitted without placement on assembly',
                                            'Submitted on non-INSDC assembly') THEN 'no_position'
        WHEN m.variant_placement_method LIKE 'Remapped%' THEN 'remapped'
        WHEN (m.variant_placement_method IS NULL OR m.variant_placement_method LIKE 'Submitted%')
             AND rm.variant_placement_id IS NULL THEN 'submitted'
        ELSE 'other'
    END AS placement_class,
    CASE WHEN EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
         THEN 'Y' ELSE 'N' END AS withdrawn,
    CASE WHEN vs.region_or_call <> 'call' THEN NULL
         WHEN EXISTS (SELECT 1 FROM dgva.linked_variant lv
                          JOIN dgva.variant_study rs ON rs.variant_accession = lv.variant_region_accession
                      WHERE lv.linked_variant_accession = vs.variant_accession
                          AND rs.study_accession = vs.study_accession) THEN 'Y'
         ELSE 'N' END AS linked_same_study,
    CASE WHEN COALESCE(st.start_pos, os.outer_start_pos, i_s.inner_start_pos) IS NULL THEN 'N' ELSE 'Y' END
        AS has_start,
    CASE WHEN REGEXP_REPLACE(p.chromosome_name, 'NOT_GIVEN', '') IS NOT NULL THEN 'Y'
         WHEN EXISTS (SELECT 1 FROM dgva.variant_placement_contig ct
                      WHERE ct.variant_placement_id = p.variant_placement_id AND ct.contig_id IS NOT NULL) THEN 'Y'
         ELSE 'N' END AS has_name
FROM dgva.variant_study vs
    LEFT JOIN dgva.variant_placement p ON p.variant_study_id = vs.variant_study_id
    LEFT JOIN dgva.genome_assembly ga ON ga.assembly_accession = p.assembly_accession
    LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_remap rm ON rm.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_start st ON st.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_outer_start os ON os.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_inner_start i_s ON i_s.variant_placement_id = p.variant_placement_id
WHERE vs.study_accession = :study
ORDER BY vs.variant_accession, ga.assembly_name
