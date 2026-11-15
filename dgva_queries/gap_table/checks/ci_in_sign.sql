-- Ticket 7, open check (gap rows 96-99): sign convention of VARIANT_PLACEMENT_CI_IN_START / _CI_IN_STOP.
-- One row per placement that has either table (845 rows), with every coordinate of the placement and its method.
-- No binds.
SELECT vs.study_accession, vs.variant_accession, vs.region_or_call, p.variant_placement_id,
    NVL(m.variant_placement_method, '(null)') AS placement_method,
    s.start_pos, os.outer_start_pos, ist.inner_start_pos, cis.ci_minus AS in_start_ci_minus, cis.ci_plus AS in_start_ci_plus,
    e.stop_pos, ie.inner_stop_pos, oe.outer_stop_pos, cie.ci_minus AS in_stop_ci_minus, cie.ci_plus AS in_stop_ci_plus
FROM dgva.variant_placement p
    JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
    LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_ci_in_start cis ON cis.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_ci_in_stop cie ON cie.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_start s ON s.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_outer_start os ON os.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_inner_start ist ON ist.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_stop e ON e.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_inner_stop ie ON ie.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_outer_stop oe ON oe.variant_placement_id = p.variant_placement_id
WHERE p.variant_placement_id IN (SELECT variant_placement_id FROM dgva.variant_placement_ci_in_start
                                 UNION SELECT variant_placement_id FROM dgva.variant_placement_ci_in_stop)
ORDER BY vs.variant_accession, p.variant_placement_id
