-- Submitted placements of one study with their CI_START/CI_STOP and CIPOS/CIEND values and their start, stop,
-- inner and outer coordinates.
-- Binds: study (one study accession).
SELECT
    vs.variant_accession,
    vs.region_or_call,
    st.start_pos, i_s.inner_start_pos, os.outer_start_pos,
    cs.ci_minus AS ci_start_minus, cs.ci_plus AS ci_start_plus,
    pl.ciposleft, pr.ciposright,
    sp.stop_pos, i_p.inner_stop_pos, op.outer_stop_pos,
    ce.ci_minus AS ci_stop_minus, ce.ci_plus AS ci_stop_plus,
    el.ciendleft, er.ciendright
FROM dgva.variant_study vs
    JOIN dgva.variant_placement p ON p.variant_study_id = vs.variant_study_id
    LEFT JOIN dgva.variant_placement_remap rm ON rm.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_ci_start cs ON cs.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_ci_stop ce ON ce.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_ciposleft pl ON pl.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_ciposright pr ON pr.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_ciendleft el ON el.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_ciendright er ON er.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_start st ON st.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_stop sp ON sp.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_inner_start i_s ON i_s.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_outer_start os ON os.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_inner_stop i_p ON i_p.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_outer_stop op ON op.variant_placement_id = p.variant_placement_id
WHERE vs.study_accession = :study
    AND rm.variant_placement_id IS NULL
    AND (cs.variant_placement_id IS NOT NULL OR pl.variant_placement_id IS NOT NULL)
