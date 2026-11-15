-- Ticket 7. Variant-level candidate gap columns (reference: candidate gap list), counted on SUBMITTED variants:
-- not withdrawn, with at least one submitted placement (ticket 2 rule). One row per study and variant type.
-- <name>_v = values (rows with a non-null value), <name>_var = variants with at least one value.
-- Tables that can hold several rows per variant are aggregated per accession first.
-- Binds: study (a study accession, or ALL). Try one study first. Reads the variant tables once: run outside busy hours.
WITH sv AS (
    SELECT vs.variant_study_id, vs.study_accession, vs.variant_accession, vs.region_or_call,
        vs.variant_load_version, vs.created_date
    FROM dgva.variant_study vs
    WHERE (:study = 'ALL' OR vs.study_accession = :study)
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
        AND EXISTS (
            SELECT 1
            FROM dgva.variant_placement p
                LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
            WHERE p.variant_study_id = vs.variant_study_id
                AND NOT EXISTS (SELECT 1 FROM dgva.variant_placement_remap rm
                                WHERE rm.variant_placement_id = p.variant_placement_id)
                AND (m.variant_placement_method IS NULL
                     OR (m.variant_placement_method LIKE 'Submitted%'
                         AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                                'Submitted on non-INSDC assembly'))))
),
call_seq AS (
    SELECT variant_call_accession AS acc, SUM(CASE WHEN variant_call_seq IS NOT NULL THEN 1 ELSE 0 END) AS n
    FROM dgva.variant_call_sequence GROUP BY variant_call_accession
),
ref_cn AS (
    SELECT variant_call_accession AS acc, COUNT(reference_copy_number) AS n
    FROM dgva.variant_call_copy_number GROUP BY variant_call_accession
),
call_sampleset AS (
    SELECT variant_call_accession AS acc, COUNT(sampleset_id) AS n
    FROM dgva.variant_call_sampleset GROUP BY variant_call_accession
),
validation AS (
    SELECT variant_call_accession AS acc, COUNT(validation_experiment_id) AS n_exp, COUNT(validation_result) AS n_res
    FROM dgva.variant_call_validation GROUP BY variant_call_accession
),
call_repeat AS (
    SELECT variant_call_accession AS acc, COUNT(repeat_count) AS n
    FROM dgva.var_call_repeat_count GROUP BY variant_call_accession
),
reg_repeat_counts AS (
    SELECT variant_region_accession AS acc, COUNT(repeat_counts) AS n
    FROM dgva.var_reg_repeat_counts GROUP BY variant_region_accession
),
reg_repeat_motif AS (
    SELECT variant_region_accession AS acc, COUNT(repeat_motif) AS n
    FROM dgva.var_reg_repeat_motif GROUP BY variant_region_accession
),
link_mutation_order AS (
    SELECT lv.variant_region_accession AS acc, COUNT(mo.mutation_order) AS n
    FROM dgva.linked_variant lv
        JOIN dgva.linked_variant_mutation_order mo ON mo.linked_variant_id = lv.linked_variant_id
    GROUP BY lv.variant_region_accession
),
call_ext_link AS (
    SELECT cl.variant_call_accession AS acc, COUNT(el.linked_object_type) AS n_type,
        COUNT(el.var_call_support_id) AS n_support_id
    FROM dgva.variant_call_ext_link cl
        JOIN dgva.external_link el ON el.external_link_id = cl.external_link_id
    GROUP BY cl.variant_call_accession
),
call_support_link AS (
    SELECT variant_call_accession AS acc, COUNT(*) AS n
    FROM dgva.var_call_support_ext_link GROUP BY variant_call_accession
),
per_variant AS (
    SELECT
        sv.study_accession,
        sv.region_or_call,
        CASE WHEN c.submitter_variant_call_type IS NOT NULL THEN 1 ELSE 0 END AS submitter_call_type,
        CASE WHEN r.submitter_variant_region_type IS NOT NULL THEN 1 ELSE 0 END AS submitter_region_type,
        CASE WHEN r.exchange_assertion_method_id IS NOT NULL THEN 1 ELSE 0 END AS exchange_assertion_method,
        CASE WHEN sv.variant_load_version > 1 THEN 1 ELSE 0 END AS load_version_above_1,
        CASE WHEN sv.created_date IS NOT NULL THEN 1 ELSE 0 END AS created_date,
        COALESCE(call_seq.n, 0) AS call_seq,
        COALESCE(ref_cn.n, 0) AS reference_copy_number,
        COALESCE(call_sampleset.n, 0) AS call_sampleset,
        COALESCE(validation.n_exp, 0) AS validation_experiment,
        COALESCE(validation.n_res, 0) AS validation_result,
        COALESCE(call_repeat.n, 0) AS call_repeat_count,
        COALESCE(reg_repeat_counts.n, 0) AS region_repeat_counts,
        COALESCE(reg_repeat_motif.n, 0) AS region_repeat_motif,
        COALESCE(link_mutation_order.n, 0) AS link_mutation_order,
        COALESCE(call_ext_link.n_type, 0) AS ext_link_object_type,
        COALESCE(call_ext_link.n_support_id, 0) AS ext_link_support_id,
        COALESCE(call_support_link.n, 0) AS call_support_link
    FROM sv
        LEFT JOIN dgva.variant_call c ON sv.region_or_call = 'call' AND c.variant_call_accession = sv.variant_accession
        LEFT JOIN dgva.variant_region r
            ON sv.region_or_call = 'region' AND r.variant_region_accession = sv.variant_accession
        LEFT JOIN call_seq ON call_seq.acc = sv.variant_accession
        LEFT JOIN ref_cn ON ref_cn.acc = sv.variant_accession
        LEFT JOIN call_sampleset ON call_sampleset.acc = sv.variant_accession
        LEFT JOIN validation ON validation.acc = sv.variant_accession
        LEFT JOIN call_repeat ON call_repeat.acc = sv.variant_accession
        LEFT JOIN reg_repeat_counts ON reg_repeat_counts.acc = sv.variant_accession
        LEFT JOIN reg_repeat_motif ON reg_repeat_motif.acc = sv.variant_accession
        LEFT JOIN link_mutation_order ON link_mutation_order.acc = sv.variant_accession
        LEFT JOIN call_ext_link ON call_ext_link.acc = sv.variant_accession
        LEFT JOIN call_support_link ON call_support_link.acc = sv.variant_accession
)
SELECT
    study_accession,
    region_or_call,
    COUNT(*) AS variants,
    SUM(submitter_call_type) AS submitter_call_type_var,
    SUM(submitter_region_type) AS submitter_region_type_var,
    SUM(exchange_assertion_method) AS exch_assertion_method_var,
    SUM(load_version_above_1) AS load_version_above_1_var,
    SUM(created_date) AS created_date_var,
    SUM(call_seq) AS call_seq_v, SUM(SIGN(call_seq)) AS call_seq_var,
    SUM(reference_copy_number) AS ref_copy_number_v, SUM(SIGN(reference_copy_number)) AS ref_copy_number_var,
    SUM(call_sampleset) AS call_sampleset_v, SUM(SIGN(call_sampleset)) AS call_sampleset_var,
    SUM(validation_experiment) AS validation_exp_v, SUM(SIGN(validation_experiment)) AS validation_exp_var,
    SUM(validation_result) AS validation_result_v, SUM(SIGN(validation_result)) AS validation_result_var,
    SUM(call_repeat_count) AS call_repeat_count_v, SUM(SIGN(call_repeat_count)) AS call_repeat_count_var,
    SUM(region_repeat_counts) AS reg_repeat_counts_v, SUM(SIGN(region_repeat_counts)) AS reg_repeat_counts_var,
    SUM(region_repeat_motif) AS reg_repeat_motif_v, SUM(SIGN(region_repeat_motif)) AS reg_repeat_motif_var,
    SUM(link_mutation_order) AS link_mutation_order_v, SUM(SIGN(link_mutation_order)) AS link_mutation_order_var,
    SUM(ext_link_object_type) AS ext_link_object_type_v, SUM(SIGN(ext_link_object_type)) AS ext_link_object_type_var,
    SUM(ext_link_support_id) AS ext_link_support_id_v, SUM(SIGN(ext_link_support_id)) AS ext_link_support_id_var,
    SUM(call_support_link) AS call_support_link_v, SUM(SIGN(call_support_link)) AS call_support_link_var
FROM per_variant
GROUP BY study_accession, region_or_call
ORDER BY study_accession, region_or_call
