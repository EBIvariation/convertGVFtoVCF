-- Every placement of the non-withdrawn variants of one study that have no submitted placement, with the variant's
-- class (same classes as ../variants_without_submitted.sql):
--   remapped_only:          remapped placements and no submitted one
--   other_placements_only:  only placements with no position
--   no_placement:           no placement at all
-- The output is the input of find_variants_in_exchange.py search.
-- Binds: study (one study accession).
WITH placements AS (
    SELECT
        vs.variant_study_id,
        vs.variant_accession,
        vs.submitter_variant_id,
        vs.region_or_call,
        p.variant_placement_id,
        p.assembly_accession,
        p.chromosome_name,
        p.chr_accession,
        m.variant_placement_method,
        r.remap_score,
        CASE
            WHEN p.variant_placement_id IS NULL THEN 'none'
            WHEN m.variant_placement_method IN ('Remap failed', 'Submitted without placement on assembly',
                                                'Submitted on non-INSDC assembly') THEN 'no_position'
            WHEN m.variant_placement_method LIKE 'Remapped%' THEN 'remapped'
            WHEN (m.variant_placement_method IS NULL OR m.variant_placement_method LIKE 'Submitted%')
                 AND r.variant_placement_id IS NULL THEN 'submitted'
            ELSE 'other'
        END AS cls
    FROM dgva.variant_study vs
        LEFT JOIN dgva.variant_placement p ON p.variant_study_id = vs.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.variant_placement_remap r ON r.variant_placement_id = p.variant_placement_id
    WHERE vs.study_accession = :study
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
),
variants AS (
    SELECT
        variant_study_id,
        CASE WHEN SUM(CASE WHEN cls = 'remapped' THEN 1 ELSE 0 END) > 0 THEN 'remapped_only'
             WHEN SUM(CASE WHEN cls = 'none' THEN 1 ELSE 0 END) > 0 THEN 'no_placement'
             ELSE 'other_placements_only' END AS variant_case
    FROM placements
    GROUP BY variant_study_id
    HAVING SUM(CASE WHEN cls = 'submitted' THEN 1 ELSE 0 END) = 0
)
SELECT
    v.variant_case,
    p.variant_accession,
    p.submitter_variant_id,
    p.region_or_call,
    p.variant_placement_id,
    p.cls AS placement_class,
    p.variant_placement_method,
    ga.assembly_name,
    p.assembly_accession,
    p.chromosome_name,
    p.chr_accession,
    st.start_pos,
    sp.stop_pos,
    os.outer_start_pos,
    op.outer_stop_pos,
    i_s.inner_start_pos,
    i_p.inner_stop_pos,
    p.remap_score,
    fc.remap_failure_code
FROM variants v
    JOIN placements p ON p.variant_study_id = v.variant_study_id
    LEFT JOIN dgva.genome_assembly ga ON ga.assembly_accession = p.assembly_accession
    LEFT JOIN dgva.variant_placement_start st ON st.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_stop sp ON sp.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_outer_start os ON os.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_outer_stop op ON op.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_inner_start i_s ON i_s.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_inner_stop i_p ON i_p.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.var_plac_remap_failure_code fc ON fc.variant_placement_id = p.variant_placement_id
ORDER BY v.variant_case, p.variant_accession, ga.assembly_name, p.variant_placement_id
