-- Per study and placement class, placements with a VARIANT_PLACEMENT_CI_START/_CI_STOP row, a
-- VARIANT_PLACEMENT_CIPOS*/CIEND* row, or both, and where both exist whether CIPOSLEFT = CI_START.CI_MINUS,
-- CIPOSRIGHT = CI_START.CI_PLUS, CIENDLEFT = CI_STOP.CI_MINUS and CIENDRIGHT = CI_STOP.CI_PLUS (missing values count
-- as 0). Also compares the intervals with the inner and outer coordinates of the same placement.
-- Placement class: submitted, remapped or no position. Submitted placement: a placement is submitted when its
-- VARIANT_PLACEMENT_METHOD is null or starts with "Submitted" (apart from "Submitted without placement on assembly" and
-- "Submitted on non-INSDC assembly") and it has no VARIANT_PLACEMENT_REMAP row.
-- Binds: none. Reads only placements that have a CI row.
WITH ci_placements AS (
    SELECT variant_placement_id FROM dgva.variant_placement_ci_start
    UNION SELECT variant_placement_id FROM dgva.variant_placement_ci_stop
    UNION SELECT variant_placement_id FROM dgva.variant_placement_ciposleft
    UNION SELECT variant_placement_id FROM dgva.variant_placement_ciposright
    UNION SELECT variant_placement_id FROM dgva.variant_placement_ciendleft
    UNION SELECT variant_placement_id FROM dgva.variant_placement_ciendright
    UNION SELECT variant_placement_id FROM dgva.variant_placement_ci_in_start
    UNION SELECT variant_placement_id FROM dgva.variant_placement_ci_in_stop
),
x AS (
    SELECT
        vs.study_accession,
        vs.region_or_call,
        CASE
            WHEN m.variant_placement_method IN ('Remap failed', 'Submitted without placement on assembly',
                                                'Submitted on non-INSDC assembly') THEN 'no_position'
            WHEN m.variant_placement_method LIKE 'Remapped%' THEN 'remapped'
            WHEN (m.variant_placement_method IS NULL OR m.variant_placement_method LIKE 'Submitted%')
                 AND rm.variant_placement_id IS NULL THEN 'submitted'
            ELSE 'other'
        END AS placement_class,
        cs.variant_placement_id AS cs_id, cs.ci_plus AS cs_plus, cs.ci_minus AS cs_minus,
        ce.variant_placement_id AS ce_id, ce.ci_plus AS ce_plus, ce.ci_minus AS ce_minus,
        pl.ciposleft, pr.ciposright, el.ciendleft, er.ciendright,
        cis.variant_placement_id AS cis_id,
        st.start_pos, sp.stop_pos, i_s.inner_start_pos, os.outer_start_pos, i_p.inner_stop_pos, op.outer_stop_pos
    FROM ci_placements c
        JOIN dgva.variant_placement p ON p.variant_placement_id = c.variant_placement_id
        JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = c.variant_placement_id
        LEFT JOIN dgva.variant_placement_remap rm ON rm.variant_placement_id = c.variant_placement_id
        LEFT JOIN dgva.variant_placement_ci_start cs ON cs.variant_placement_id = c.variant_placement_id
        LEFT JOIN dgva.variant_placement_ci_stop ce ON ce.variant_placement_id = c.variant_placement_id
        LEFT JOIN dgva.variant_placement_ciposleft pl ON pl.variant_placement_id = c.variant_placement_id
        LEFT JOIN dgva.variant_placement_ciposright pr ON pr.variant_placement_id = c.variant_placement_id
        LEFT JOIN dgva.variant_placement_ciendleft el ON el.variant_placement_id = c.variant_placement_id
        LEFT JOIN dgva.variant_placement_ciendright er ON er.variant_placement_id = c.variant_placement_id
        LEFT JOIN dgva.variant_placement_ci_in_start cis ON cis.variant_placement_id = c.variant_placement_id
        LEFT JOIN dgva.variant_placement_start st ON st.variant_placement_id = c.variant_placement_id
        LEFT JOIN dgva.variant_placement_stop sp ON sp.variant_placement_id = c.variant_placement_id
        LEFT JOIN dgva.variant_placement_inner_start i_s ON i_s.variant_placement_id = c.variant_placement_id
        LEFT JOIN dgva.variant_placement_outer_start os ON os.variant_placement_id = c.variant_placement_id
        LEFT JOIN dgva.variant_placement_inner_stop i_p ON i_p.variant_placement_id = c.variant_placement_id
        LEFT JOIN dgva.variant_placement_outer_stop op ON op.variant_placement_id = c.variant_placement_id
)
SELECT
    study_accession,
    region_or_call,
    placement_class,
    COUNT(*) AS placements_with_any_ci,
    COUNT(cs_id) AS with_ci_start,
    COUNT(ce_id) AS with_ci_stop,
    SUM(CASE WHEN cs_id IS NOT NULL AND cs_plus = 0 AND cs_minus = 0 THEN 1 ELSE 0 END) AS ci_start_all_zero,
    SUM(CASE WHEN ciposleft IS NOT NULL OR ciposright IS NOT NULL THEN 1 ELSE 0 END) AS with_cipos,
    SUM(CASE WHEN ciendleft IS NOT NULL OR ciendright IS NOT NULL THEN 1 ELSE 0 END) AS with_ciend,
    COUNT(cis_id) AS with_ci_in_start,
    SUM(CASE WHEN cs_id IS NOT NULL AND (ciposleft IS NOT NULL OR ciposright IS NOT NULL) THEN 1 ELSE 0 END)
        AS both_start_sets,
    SUM(CASE WHEN cs_id IS NOT NULL AND (ciposleft IS NOT NULL OR ciposright IS NOT NULL)
              AND COALESCE(ciposleft, 0) = cs_minus AND COALESCE(ciposright, 0) = cs_plus THEN 1 ELSE 0 END)
        AS both_start_agree,
    SUM(CASE WHEN ce_id IS NOT NULL AND (ciendleft IS NOT NULL OR ciendright IS NOT NULL) THEN 1 ELSE 0 END)
        AS both_stop_sets,
    SUM(CASE WHEN ce_id IS NOT NULL AND (ciendleft IS NOT NULL OR ciendright IS NOT NULL)
              AND COALESCE(ciendleft, 0) = ce_minus AND COALESCE(ciendright, 0) = ce_plus THEN 1 ELSE 0 END)
        AS both_stop_agree,
    SUM(CASE WHEN cs_id IS NOT NULL AND start_pos IS NOT NULL THEN 1 ELSE 0 END) AS ci_start_with_start,
    SUM(CASE WHEN cs_id IS NOT NULL AND (inner_start_pos IS NOT NULL OR outer_start_pos IS NOT NULL)
             THEN 1 ELSE 0 END) AS ci_start_with_inner_or_outer,
    SUM(CASE WHEN cs_id IS NOT NULL AND start_pos IS NOT NULL
              AND (inner_start_pos IS NOT NULL OR outer_start_pos IS NOT NULL)
              AND COALESCE(inner_start_pos, start_pos + cs_plus) = start_pos + cs_plus
              AND COALESCE(outer_start_pos, start_pos - cs_minus) = start_pos - cs_minus THEN 1 ELSE 0 END)
        AS ci_start_matches_inner_outer,
    SUM(CASE WHEN ce_id IS NOT NULL AND stop_pos IS NOT NULL
              AND (inner_stop_pos IS NOT NULL OR outer_stop_pos IS NOT NULL)
              AND COALESCE(inner_stop_pos, stop_pos - ce_minus) = stop_pos - ce_minus
              AND COALESCE(outer_stop_pos, stop_pos + ce_plus) = stop_pos + ce_plus THEN 1 ELSE 0 END)
        AS ci_stop_matches_inner_outer,
    SUM(CASE WHEN cs_id IS NULL AND ce_id IS NULL AND start_pos IS NOT NULL
              AND (ciposleft IS NOT NULL OR ciposright IS NOT NULL) THEN 1 ELSE 0 END) AS cipos_only
FROM x
GROUP BY study_accession, region_or_call, placement_class
ORDER BY study_accession, region_or_call, placement_class
