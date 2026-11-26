-- Per study, assembly and variant kind: submitted placements of non-withdrawn variants, flagged with the conditions
-- under which archive_gvf/get_gvf_from_archive.pl and table_queries.sql (dgva-pipeline) write a placement:
--   W2  the variant is in VARIANT_REGION or VARIANT_CALL, its type maps to a SEQUENCE_ONTOLOGY_CV term, a region has a
--       VAR_REG_ASSERT_METHOD row, and the assembly is in GENOME_ASSEMBLY;
--   W1  a call is linked (LINKED_VARIANT) from a region of the same study;
--   W3  the placement has a start (start, outer start or inner start) and a chromosome or contig name.
-- Submitted placement: a placement is submitted when its VARIANT_PLACEMENT_METHOD is null or starts with "Submitted"
-- (apart from "Submitted without placement on assembly" and "Submitted on non-INSDC assembly") and it has no
-- VARIANT_PLACEMENT_REMAP row.
-- Binds: study (a study accession, or ALL). With ALL, reads every placement.
WITH linked_same_study AS (
    SELECT DISTINCT lv.linked_variant_accession AS call_accession, rs.study_accession
    FROM dgva.linked_variant lv
        JOIN dgva.variant_region vr ON vr.variant_region_accession = lv.variant_region_accession
        JOIN dgva.variant_study rs ON rs.variant_accession = lv.variant_region_accession
    WHERE (:study = 'ALL' OR rs.study_accession = :study)
),
linked_any AS (
    SELECT DISTINCT lv.linked_variant_accession AS call_accession
    FROM dgva.linked_variant lv
),
v AS (
    SELECT
        vs.variant_study_id,
        vs.study_accession,
        vs.region_or_call,
        CASE WHEN COALESCE(r.variant_region_accession, c.variant_call_accession) IS NULL THEN 'N' ELSE 'Y' END
            AS in_type_table,
        CASE
            WHEN vs.region_or_call = 'region' AND EXISTS (
                SELECT 1 FROM dgva.variant_region_type_cv t
                    JOIN dgva.sequence_ontology_cv so ON so.sequence_ontology_term = t.sequence_ontology_term
                WHERE t.exchange_variant_region_type = r.exchange_variant_region_type) THEN 'Y'
            WHEN vs.region_or_call = 'call' AND EXISTS (
                SELECT 1 FROM dgva.variant_call_type_cv t
                    JOIN dgva.sequence_ontology_cv so ON so.sequence_ontology_term = t.sequence_ontology_term
                WHERE t.exchange_variant_call_type = c.exchange_variant_call_type) THEN 'Y'
            ELSE 'N'
        END AS has_so,
        CASE WHEN vs.region_or_call <> 'region' THEN 'Y'
             WHEN EXISTS (SELECT 1 FROM dgva.var_reg_assert_method a
                          WHERE a.var_reg_assert_method_id = r.var_reg_assert_method_id) THEN 'Y'
             ELSE 'N'
        END AS has_assertion,
        CASE WHEN la.call_accession IS NULL THEN 'N' ELSE 'Y' END AS linked_any,
        CASE WHEN ls.call_accession IS NULL THEN 'N' ELSE 'Y' END AS linked_same_study
    FROM dgva.variant_study vs
        LEFT JOIN dgva.variant_region r
            ON vs.region_or_call = 'region' AND r.variant_region_accession = vs.variant_accession
        LEFT JOIN dgva.variant_call c
            ON vs.region_or_call = 'call' AND c.variant_call_accession = vs.variant_accession
        LEFT JOIN linked_any la
            ON vs.region_or_call = 'call' AND la.call_accession = vs.variant_accession
        LEFT JOIN linked_same_study ls
            ON vs.region_or_call = 'call' AND ls.call_accession = vs.variant_accession
            AND ls.study_accession = vs.study_accession
    WHERE (:study = 'ALL' OR vs.study_accession = :study)
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
),
p AS (
    SELECT
        p.variant_placement_id,
        p.variant_study_id,
        p.assembly_accession,
        CASE WHEN EXISTS (SELECT 1 FROM dgva.genome_assembly ga WHERE ga.assembly_accession = p.assembly_accession)
             THEN 'Y' ELSE 'N' END AS has_assembly,
        CASE WHEN COALESCE(st.start_pos, os.outer_start_pos, i_s.inner_start_pos) IS NULL THEN 'N' ELSE 'Y' END
            AS has_start,
        CASE WHEN REGEXP_REPLACE(p.chromosome_name, 'NOT_GIVEN', '') IS NOT NULL THEN 'Y'
             WHEN EXISTS (SELECT 1 FROM dgva.variant_placement_contig ct
                          WHERE ct.variant_placement_id = p.variant_placement_id AND ct.contig_id IS NOT NULL) THEN 'Y'
             ELSE 'N'
        END AS has_name
    FROM dgva.variant_placement p
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.variant_placement_remap rm ON rm.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.variant_placement_start st ON st.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.variant_placement_outer_start os ON os.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.variant_placement_inner_start i_s ON i_s.variant_placement_id = p.variant_placement_id
    WHERE rm.variant_placement_id IS NULL
        AND (m.variant_placement_method IS NULL
             OR (m.variant_placement_method LIKE 'Submitted%'
                 AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                        'Submitted on non-INSDC assembly')))
),
pv AS (
    SELECT
        v.*,
        p.variant_placement_id,
        p.assembly_accession,
        p.has_assembly,
        p.has_start,
        p.has_name,
        CASE WHEN v.in_type_table = 'Y' AND v.has_so = 'Y' AND v.has_assertion = 'Y'
                  AND (v.region_or_call = 'region' OR v.linked_same_study = 'Y')
                  AND p.has_assembly = 'Y' AND p.has_start = 'Y' AND p.has_name = 'Y'
             THEN 'Y' ELSE 'N' END AS written
    FROM v
        JOIN p ON p.variant_study_id = v.variant_study_id
)
SELECT
    pv.study_accession,
    ga.assembly_name,
    pv.assembly_accession,
    pv.region_or_call,
    COUNT(DISTINCT pv.variant_study_id) AS variants,
    COUNT(DISTINCT CASE WHEN pv.written = 'Y' THEN pv.variant_study_id END) AS variants_written,
    COUNT(DISTINCT pv.variant_study_id)
        - COUNT(DISTINCT CASE WHEN pv.written = 'Y' THEN pv.variant_study_id END) AS variants_dropped,
    COUNT(DISTINCT CASE WHEN pv.region_or_call = 'call' AND pv.linked_any = 'N'
                        THEN pv.variant_study_id END) AS w1_calls_no_region,
    COUNT(DISTINCT CASE WHEN pv.region_or_call = 'call' AND pv.linked_any = 'Y' AND pv.linked_same_study = 'N'
                        THEN pv.variant_study_id END) AS w1_calls_other_study_region_only,
    COUNT(DISTINCT CASE WHEN pv.in_type_table = 'N' THEN pv.variant_study_id END) AS w2_not_in_type_table,
    COUNT(DISTINCT CASE WHEN pv.in_type_table = 'Y' AND pv.has_so = 'N' THEN pv.variant_study_id END)
        AS w2_type_without_so,
    COUNT(DISTINCT CASE WHEN pv.has_assertion = 'N' THEN pv.variant_study_id END) AS w2_regions_no_assertion,
    COUNT(*) AS placements,
    SUM(CASE WHEN pv.written = 'Y' THEN 1 ELSE 0 END) AS placements_written,
    SUM(CASE WHEN pv.has_assembly = 'N' THEN 1 ELSE 0 END) AS w2_placements_assembly_unknown,
    SUM(CASE WHEN pv.has_start = 'N' THEN 1 ELSE 0 END) AS w3_placements_no_start,
    SUM(CASE WHEN pv.has_start = 'Y' AND pv.has_name = 'N' THEN 1 ELSE 0 END) AS w3_placements_start_but_no_name
FROM pv
    LEFT JOIN dgva.genome_assembly ga ON ga.assembly_accession = pv.assembly_accession
GROUP BY pv.study_accession, ga.assembly_name, pv.assembly_accession, pv.region_or_call
ORDER BY pv.study_accession, ga.assembly_name, pv.region_or_call
