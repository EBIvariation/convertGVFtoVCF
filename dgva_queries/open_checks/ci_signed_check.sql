-- Ticket 8.1. For CI_START/CI_STOP, test the signed reading seen in the data (CI_MINUS stored negative, as VCF CIPOS):
-- outer_start = start + CI_START.CI_MINUS, inner_start = start + CI_START.CI_PLUS,
-- inner_stop = stop + CI_STOP.CI_MINUS, outer_stop = stop + CI_STOP.CI_PLUS.
-- If every placement matches, the inner/outer coordinates already carry the same intervals.
-- Binds: study (one study accession). One row per variant type.
SELECT
    vs.region_or_call,
    COUNT(*) AS placements_with_ci_start,
    SUM(CASE WHEN cs.ci_minus > 0 OR ce.ci_minus > 0 THEN 1 ELSE 0 END) AS ci_minus_positive,
    SUM(CASE WHEN DECODE(os.outer_start_pos, st.start_pos + cs.ci_minus, 1, 0) = 1
              AND DECODE(i_s.inner_start_pos, st.start_pos + cs.ci_plus, 1, 0) = 1 THEN 1 ELSE 0 END) AS start_matches,
    SUM(CASE WHEN DECODE(i_p.inner_stop_pos, sp.stop_pos + ce.ci_minus, 1, 0) = 1
              AND DECODE(op.outer_stop_pos, sp.stop_pos + ce.ci_plus, 1, 0) = 1 THEN 1 ELSE 0 END) AS stop_matches
FROM dgva.variant_study vs
    JOIN dgva.variant_placement p ON p.variant_study_id = vs.variant_study_id
    JOIN dgva.variant_placement_ci_start cs ON cs.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_ci_stop ce ON ce.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_start st ON st.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_stop sp ON sp.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_inner_start i_s ON i_s.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_outer_start os ON os.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_inner_stop i_p ON i_p.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_outer_stop op ON op.variant_placement_id = p.variant_placement_id
WHERE vs.study_accession = :study
GROUP BY vs.region_or_call
