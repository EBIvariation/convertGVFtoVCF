-- Per table, column, study and variant type: submitted placements of non-withdrawn variants with a value in the
-- placement-level column, and variants with a value on at least one such placement.
-- Submitted placement: a placement is submitted when its VARIANT_PLACEMENT_METHOD is null or starts with "Submitted"
-- (apart from "Submitted without placement on assembly" and "Submitted on non-INSDC assembly") and it has no
-- VARIANT_PLACEMENT_REMAP row.
-- Each branch runs only when :table_name is its table or 'ALL'. Example: --bind table_name=VAR_PLAC_PLACEMENT_RANK
-- --bind study=ALL.
-- Binds: table_name (a table below, or ALL), study (a study accession, or ALL).
SELECT * FROM (
    SELECT 'VARIANT_PLACEMENT' AS table_name, 'CHR_ACCESSION' AS column_name, vs.study_accession, vs.region_or_call,
        COUNT(*) AS placements_with_value, COUNT(DISTINCT vs.variant_study_id) AS variants_with_value
    FROM dgva.variant_placement p
        JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
    WHERE :table_name IN ('ALL', 'VARIANT_PLACEMENT')
        AND (:study = 'ALL' OR vs.study_accession = :study)
        AND p.chr_accession IS NOT NULL
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_placement_remap rm WHERE rm.variant_placement_id = p.variant_placement_id)
        AND (m.variant_placement_method IS NULL
             OR (m.variant_placement_method LIKE 'Submitted%'
                 AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                        'Submitted on non-INSDC assembly')))
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
    GROUP BY vs.study_accession, vs.region_or_call
    UNION ALL
    SELECT 'VARIANT_PLACEMENT_ALT_STATUS' AS table_name, 'ALT_STATUS' AS column_name, vs.study_accession, vs.region_or_call,
        COUNT(*) AS placements_with_value, COUNT(DISTINCT vs.variant_study_id) AS variants_with_value
    FROM dgva.variant_placement_alt_status t
        JOIN dgva.variant_placement p ON p.variant_placement_id = t.variant_placement_id
        JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
    WHERE :table_name IN ('ALL', 'VARIANT_PLACEMENT_ALT_STATUS')
        AND (:study = 'ALL' OR vs.study_accession = :study)
        AND t.alt_status IS NOT NULL
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_placement_remap rm WHERE rm.variant_placement_id = p.variant_placement_id)
        AND (m.variant_placement_method IS NULL
             OR (m.variant_placement_method LIKE 'Submitted%'
                 AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                        'Submitted on non-INSDC assembly')))
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
    GROUP BY vs.study_accession, vs.region_or_call
    UNION ALL
    SELECT 'VARIANT_PLACEMENT_ASSEM_UNIT' AS table_name, 'ASSEMBLY_UNIT' AS column_name, vs.study_accession, vs.region_or_call,
        COUNT(*) AS placements_with_value, COUNT(DISTINCT vs.variant_study_id) AS variants_with_value
    FROM dgva.variant_placement_assem_unit t
        JOIN dgva.variant_placement p ON p.variant_placement_id = t.variant_placement_id
        JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
    WHERE :table_name IN ('ALL', 'VARIANT_PLACEMENT_ASSEM_UNIT')
        AND (:study = 'ALL' OR vs.study_accession = :study)
        AND t.assembly_unit IS NOT NULL
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_placement_remap rm WHERE rm.variant_placement_id = p.variant_placement_id)
        AND (m.variant_placement_method IS NULL
             OR (m.variant_placement_method LIKE 'Submitted%'
                 AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                        'Submitted on non-INSDC assembly')))
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
    GROUP BY vs.study_accession, vs.region_or_call
    UNION ALL
    SELECT 'VARIANT_PLACEMENT_CIPOSLEFT' AS table_name, 'CIPOSLEFT' AS column_name, vs.study_accession, vs.region_or_call,
        COUNT(*) AS placements_with_value, COUNT(DISTINCT vs.variant_study_id) AS variants_with_value
    FROM dgva.variant_placement_ciposleft t
        JOIN dgva.variant_placement p ON p.variant_placement_id = t.variant_placement_id
        JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
    WHERE :table_name IN ('ALL', 'VARIANT_PLACEMENT_CIPOSLEFT')
        AND (:study = 'ALL' OR vs.study_accession = :study)
        AND t.ciposleft IS NOT NULL
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_placement_remap rm WHERE rm.variant_placement_id = p.variant_placement_id)
        AND (m.variant_placement_method IS NULL
             OR (m.variant_placement_method LIKE 'Submitted%'
                 AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                        'Submitted on non-INSDC assembly')))
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
    GROUP BY vs.study_accession, vs.region_or_call
    UNION ALL
    SELECT 'VARIANT_PLACEMENT_CIPOSRIGHT' AS table_name, 'CIPOSRIGHT' AS column_name, vs.study_accession, vs.region_or_call,
        COUNT(*) AS placements_with_value, COUNT(DISTINCT vs.variant_study_id) AS variants_with_value
    FROM dgva.variant_placement_ciposright t
        JOIN dgva.variant_placement p ON p.variant_placement_id = t.variant_placement_id
        JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
    WHERE :table_name IN ('ALL', 'VARIANT_PLACEMENT_CIPOSRIGHT')
        AND (:study = 'ALL' OR vs.study_accession = :study)
        AND t.ciposright IS NOT NULL
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_placement_remap rm WHERE rm.variant_placement_id = p.variant_placement_id)
        AND (m.variant_placement_method IS NULL
             OR (m.variant_placement_method LIKE 'Submitted%'
                 AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                        'Submitted on non-INSDC assembly')))
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
    GROUP BY vs.study_accession, vs.region_or_call
    UNION ALL
    SELECT 'VARIANT_PLACEMENT_CIENDLEFT' AS table_name, 'CIENDLEFT' AS column_name, vs.study_accession, vs.region_or_call,
        COUNT(*) AS placements_with_value, COUNT(DISTINCT vs.variant_study_id) AS variants_with_value
    FROM dgva.variant_placement_ciendleft t
        JOIN dgva.variant_placement p ON p.variant_placement_id = t.variant_placement_id
        JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
    WHERE :table_name IN ('ALL', 'VARIANT_PLACEMENT_CIENDLEFT')
        AND (:study = 'ALL' OR vs.study_accession = :study)
        AND t.ciendleft IS NOT NULL
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_placement_remap rm WHERE rm.variant_placement_id = p.variant_placement_id)
        AND (m.variant_placement_method IS NULL
             OR (m.variant_placement_method LIKE 'Submitted%'
                 AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                        'Submitted on non-INSDC assembly')))
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
    GROUP BY vs.study_accession, vs.region_or_call
    UNION ALL
    SELECT 'VARIANT_PLACEMENT_CIENDRIGHT' AS table_name, 'CIENDRIGHT' AS column_name, vs.study_accession, vs.region_or_call,
        COUNT(*) AS placements_with_value, COUNT(DISTINCT vs.variant_study_id) AS variants_with_value
    FROM dgva.variant_placement_ciendright t
        JOIN dgva.variant_placement p ON p.variant_placement_id = t.variant_placement_id
        JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
    WHERE :table_name IN ('ALL', 'VARIANT_PLACEMENT_CIENDRIGHT')
        AND (:study = 'ALL' OR vs.study_accession = :study)
        AND t.ciendright IS NOT NULL
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_placement_remap rm WHERE rm.variant_placement_id = p.variant_placement_id)
        AND (m.variant_placement_method IS NULL
             OR (m.variant_placement_method LIKE 'Submitted%'
                 AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                        'Submitted on non-INSDC assembly')))
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
    GROUP BY vs.study_accession, vs.region_or_call
    UNION ALL
    SELECT 'VARIANT_PLACEMENT_CI_IN_START' AS table_name, 'CI_PLUS' AS column_name, vs.study_accession, vs.region_or_call,
        COUNT(*) AS placements_with_value, COUNT(DISTINCT vs.variant_study_id) AS variants_with_value
    FROM dgva.variant_placement_ci_in_start t
        JOIN dgva.variant_placement p ON p.variant_placement_id = t.variant_placement_id
        JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
    WHERE :table_name IN ('ALL', 'VARIANT_PLACEMENT_CI_IN_START')
        AND (:study = 'ALL' OR vs.study_accession = :study)
        AND t.ci_plus IS NOT NULL
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_placement_remap rm WHERE rm.variant_placement_id = p.variant_placement_id)
        AND (m.variant_placement_method IS NULL
             OR (m.variant_placement_method LIKE 'Submitted%'
                 AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                        'Submitted on non-INSDC assembly')))
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
    GROUP BY vs.study_accession, vs.region_or_call
    UNION ALL
    SELECT 'VARIANT_PLACEMENT_CI_IN_STOP' AS table_name, 'CI_PLUS' AS column_name, vs.study_accession, vs.region_or_call,
        COUNT(*) AS placements_with_value, COUNT(DISTINCT vs.variant_study_id) AS variants_with_value
    FROM dgva.variant_placement_ci_in_stop t
        JOIN dgva.variant_placement p ON p.variant_placement_id = t.variant_placement_id
        JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
    WHERE :table_name IN ('ALL', 'VARIANT_PLACEMENT_CI_IN_STOP')
        AND (:study = 'ALL' OR vs.study_accession = :study)
        AND t.ci_plus IS NOT NULL
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_placement_remap rm WHERE rm.variant_placement_id = p.variant_placement_id)
        AND (m.variant_placement_method IS NULL
             OR (m.variant_placement_method LIKE 'Submitted%'
                 AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                        'Submitted on non-INSDC assembly')))
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
    GROUP BY vs.study_accession, vs.region_or_call
    UNION ALL
    SELECT 'VARIANT_PLACEMENT_CI_START' AS table_name, 'CI_PLUS' AS column_name, vs.study_accession, vs.region_or_call,
        COUNT(*) AS placements_with_value, COUNT(DISTINCT vs.variant_study_id) AS variants_with_value
    FROM dgva.variant_placement_ci_start t
        JOIN dgva.variant_placement p ON p.variant_placement_id = t.variant_placement_id
        JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
    WHERE :table_name IN ('ALL', 'VARIANT_PLACEMENT_CI_START')
        AND (:study = 'ALL' OR vs.study_accession = :study)
        AND t.ci_plus IS NOT NULL
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_placement_remap rm WHERE rm.variant_placement_id = p.variant_placement_id)
        AND (m.variant_placement_method IS NULL
             OR (m.variant_placement_method LIKE 'Submitted%'
                 AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                        'Submitted on non-INSDC assembly')))
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
    GROUP BY vs.study_accession, vs.region_or_call
    UNION ALL
    SELECT 'VARIANT_PLACEMENT_CI_STOP' AS table_name, 'CI_PLUS' AS column_name, vs.study_accession, vs.region_or_call,
        COUNT(*) AS placements_with_value, COUNT(DISTINCT vs.variant_study_id) AS variants_with_value
    FROM dgva.variant_placement_ci_stop t
        JOIN dgva.variant_placement p ON p.variant_placement_id = t.variant_placement_id
        JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
    WHERE :table_name IN ('ALL', 'VARIANT_PLACEMENT_CI_STOP')
        AND (:study = 'ALL' OR vs.study_accession = :study)
        AND t.ci_plus IS NOT NULL
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_placement_remap rm WHERE rm.variant_placement_id = p.variant_placement_id)
        AND (m.variant_placement_method IS NULL
             OR (m.variant_placement_method LIKE 'Submitted%'
                 AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                        'Submitted on non-INSDC assembly')))
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
    GROUP BY vs.study_accession, vs.region_or_call
    UNION ALL
    SELECT 'VARIANT_PLACEMENT_CYTOBAND' AS table_name, 'CYTOBAND' AS column_name, vs.study_accession, vs.region_or_call,
        COUNT(*) AS placements_with_value, COUNT(DISTINCT vs.variant_study_id) AS variants_with_value
    FROM dgva.variant_placement_cytoband t
        JOIN dgva.variant_placement p ON p.variant_placement_id = t.variant_placement_id
        JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
    WHERE :table_name IN ('ALL', 'VARIANT_PLACEMENT_CYTOBAND')
        AND (:study = 'ALL' OR vs.study_accession = :study)
        AND t.cytoband IS NOT NULL
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_placement_remap rm WHERE rm.variant_placement_id = p.variant_placement_id)
        AND (m.variant_placement_method IS NULL
             OR (m.variant_placement_method LIKE 'Submitted%'
                 AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                        'Submitted on non-INSDC assembly')))
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
    GROUP BY vs.study_accession, vs.region_or_call
    UNION ALL
    SELECT 'VARIANT_PLACEMENT_CYTOBAND' AS table_name, 'CYTOBAND_TAX_ID' AS column_name, vs.study_accession, vs.region_or_call,
        COUNT(*) AS placements_with_value, COUNT(DISTINCT vs.variant_study_id) AS variants_with_value
    FROM dgva.variant_placement_cytoband t
        JOIN dgva.variant_placement p ON p.variant_placement_id = t.variant_placement_id
        JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
    WHERE :table_name IN ('ALL', 'VARIANT_PLACEMENT_CYTOBAND')
        AND (:study = 'ALL' OR vs.study_accession = :study)
        AND t.cytoband_tax_id IS NOT NULL
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_placement_remap rm WHERE rm.variant_placement_id = p.variant_placement_id)
        AND (m.variant_placement_method IS NULL
             OR (m.variant_placement_method LIKE 'Submitted%'
                 AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                        'Submitted on non-INSDC assembly')))
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
    GROUP BY vs.study_accession, vs.region_or_call
    UNION ALL
    SELECT 'VARIANT_PLACEMENT_HGVS' AS table_name, 'HGVS_NAME' AS column_name, vs.study_accession, vs.region_or_call,
        COUNT(*) AS placements_with_value, COUNT(DISTINCT vs.variant_study_id) AS variants_with_value
    FROM dgva.variant_placement_hgvs t
        JOIN dgva.variant_placement p ON p.variant_placement_id = t.variant_placement_id
        JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
    WHERE :table_name IN ('ALL', 'VARIANT_PLACEMENT_HGVS')
        AND (:study = 'ALL' OR vs.study_accession = :study)
        AND t.hgvs_name IS NOT NULL
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_placement_remap rm WHERE rm.variant_placement_id = p.variant_placement_id)
        AND (m.variant_placement_method IS NULL
             OR (m.variant_placement_method LIKE 'Submitted%'
                 AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                        'Submitted on non-INSDC assembly')))
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
    GROUP BY vs.study_accession, vs.region_or_call
    UNION ALL
    SELECT 'VAR_PLAC_MUTATION_ORDER' AS table_name, 'MUTATION_ORDER' AS column_name, vs.study_accession, vs.region_or_call,
        COUNT(*) AS placements_with_value, COUNT(DISTINCT vs.variant_study_id) AS variants_with_value
    FROM dgva.var_plac_mutation_order t
        JOIN dgva.variant_placement p ON p.variant_placement_id = t.variant_placement_id
        JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
    WHERE :table_name IN ('ALL', 'VAR_PLAC_MUTATION_ORDER')
        AND (:study = 'ALL' OR vs.study_accession = :study)
        AND t.mutation_order IS NOT NULL
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_placement_remap rm WHERE rm.variant_placement_id = p.variant_placement_id)
        AND (m.variant_placement_method IS NULL
             OR (m.variant_placement_method LIKE 'Submitted%'
                 AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                        'Submitted on non-INSDC assembly')))
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
    GROUP BY vs.study_accession, vs.region_or_call
    UNION ALL
    SELECT 'VAR_PLAC_PLACEMENT_RANK' AS table_name, 'PLACEMENT_RANK' AS column_name, vs.study_accession, vs.region_or_call,
        COUNT(*) AS placements_with_value, COUNT(DISTINCT vs.variant_study_id) AS variants_with_value
    FROM dgva.var_plac_placement_rank t
        JOIN dgva.variant_placement p ON p.variant_placement_id = t.variant_placement_id
        JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
    WHERE :table_name IN ('ALL', 'VAR_PLAC_PLACEMENT_RANK')
        AND (:study = 'ALL' OR vs.study_accession = :study)
        AND t.placement_rank IS NOT NULL
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_placement_remap rm WHERE rm.variant_placement_id = p.variant_placement_id)
        AND (m.variant_placement_method IS NULL
             OR (m.variant_placement_method LIKE 'Submitted%'
                 AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                        'Submitted on non-INSDC assembly')))
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
    GROUP BY vs.study_accession, vs.region_or_call
    UNION ALL
    SELECT 'VAR_PLAC_PLACS_PER_ASSEM' AS table_name, 'PLACEMENTS_PER_ASSEMBLY' AS column_name, vs.study_accession, vs.region_or_call,
        COUNT(*) AS placements_with_value, COUNT(DISTINCT vs.variant_study_id) AS variants_with_value
    FROM dgva.var_plac_placs_per_assem t
        JOIN dgva.variant_placement p ON p.variant_placement_id = t.variant_placement_id
        JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
    WHERE :table_name IN ('ALL', 'VAR_PLAC_PLACS_PER_ASSEM')
        AND (:study = 'ALL' OR vs.study_accession = :study)
        AND t.placements_per_assembly IS NOT NULL
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_placement_remap rm WHERE rm.variant_placement_id = p.variant_placement_id)
        AND (m.variant_placement_method IS NULL
             OR (m.variant_placement_method LIKE 'Submitted%'
                 AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                        'Submitted on non-INSDC assembly')))
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
    GROUP BY vs.study_accession, vs.region_or_call
    UNION ALL
    SELECT 'VAR_PLAC_REMAP_BEST_CLUST' AS table_name, 'REMAP_BEST_WITHIN_CLUSTER' AS column_name, vs.study_accession, vs.region_or_call,
        COUNT(*) AS placements_with_value, COUNT(DISTINCT vs.variant_study_id) AS variants_with_value
    FROM dgva.var_plac_remap_best_clust t
        JOIN dgva.variant_placement p ON p.variant_placement_id = t.variant_placement_id
        JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
    WHERE :table_name IN ('ALL', 'VAR_PLAC_REMAP_BEST_CLUST')
        AND (:study = 'ALL' OR vs.study_accession = :study)
        AND t.remap_best_within_cluster IS NOT NULL
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_placement_remap rm WHERE rm.variant_placement_id = p.variant_placement_id)
        AND (m.variant_placement_method IS NULL
             OR (m.variant_placement_method LIKE 'Submitted%'
                 AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                        'Submitted on non-INSDC assembly')))
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
    GROUP BY vs.study_accession, vs.region_or_call
    UNION ALL
    SELECT 'VAR_PLAC_REMAP_DIFF_CHR' AS table_name, 'REMAP_DIFF_CHR' AS column_name, vs.study_accession, vs.region_or_call,
        COUNT(*) AS placements_with_value, COUNT(DISTINCT vs.variant_study_id) AS variants_with_value
    FROM dgva.var_plac_remap_diff_chr t
        JOIN dgva.variant_placement p ON p.variant_placement_id = t.variant_placement_id
        JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
    WHERE :table_name IN ('ALL', 'VAR_PLAC_REMAP_DIFF_CHR')
        AND (:study = 'ALL' OR vs.study_accession = :study)
        AND t.remap_diff_chr IS NOT NULL
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_placement_remap rm WHERE rm.variant_placement_id = p.variant_placement_id)
        AND (m.variant_placement_method IS NULL
             OR (m.variant_placement_method LIKE 'Submitted%'
                 AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                        'Submitted on non-INSDC assembly')))
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
    GROUP BY vs.study_accession, vs.region_or_call
    UNION ALL
    SELECT 'VAR_PLAC_REMAP_FAILURE_CODE' AS table_name, 'REMAP_FAILURE_CODE' AS column_name, vs.study_accession, vs.region_or_call,
        COUNT(*) AS placements_with_value, COUNT(DISTINCT vs.variant_study_id) AS variants_with_value
    FROM dgva.var_plac_remap_failure_code t
        JOIN dgva.variant_placement p ON p.variant_placement_id = t.variant_placement_id
        JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
    WHERE :table_name IN ('ALL', 'VAR_PLAC_REMAP_FAILURE_CODE')
        AND (:study = 'ALL' OR vs.study_accession = :study)
        AND t.remap_failure_code IS NOT NULL
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_placement_remap rm WHERE rm.variant_placement_id = p.variant_placement_id)
        AND (m.variant_placement_method IS NULL
             OR (m.variant_placement_method LIKE 'Submitted%'
                 AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                        'Submitted on non-INSDC assembly')))
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
    GROUP BY vs.study_accession, vs.region_or_call
)
ORDER BY table_name, column_name, study_accession, region_or_call
