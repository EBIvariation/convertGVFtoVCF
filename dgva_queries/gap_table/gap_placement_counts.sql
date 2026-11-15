-- Ticket 7. Placement-level candidate gap columns (reference: candidate gap list), counted on SUBMITTED placements
-- (ticket 2 rule) of non-withdrawn variants. One row per study and variant type. For each column:
--   <alias>_v = placements with a value, <alias>_var = variants with a value on at least one submitted placement.
-- A table that can hold several rows per placement is aggregated first (COUNT of non-null values per placement).
-- CI_IN_*/CI_* are counted on CI_PLUS (the loader writes both values or 0). The column list follows the candidate gap list.
-- Binds: none. SUPERSEDED by gap_placement_table_counts.sql: this single scan did not finish in 10 hours
-- (stopped on 2026-11-10). Kept for reference only; do not run.
WITH sp AS (
    SELECT p.variant_placement_id, p.chr_accession, vs.variant_study_id, vs.study_accession, vs.region_or_call
    FROM dgva.variant_placement p
        JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.variant_placement_remap rm ON rm.variant_placement_id = p.variant_placement_id
    WHERE rm.variant_placement_id IS NULL
        AND (m.variant_placement_method IS NULL
             OR (m.variant_placement_method LIKE 'Submitted%'
                 AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                        'Submitted on non-INSDC assembly')))
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
),
per_placement AS (
    SELECT
        sp.study_accession, sp.region_or_call, sp.variant_study_id,
        CASE WHEN sp.chr_accession IS NOT NULL THEN 1 ELSE 0 END AS chr_accession,
        COALESCE(x_alt_status.n, 0) AS alt_status,
        COALESCE(x_assembly_unit.n, 0) AS assembly_unit,
        COALESCE(x_ciposleft.n, 0) AS ciposleft,
        COALESCE(x_ciposright.n, 0) AS ciposright,
        COALESCE(x_ciendleft.n, 0) AS ciendleft,
        COALESCE(x_ciendright.n, 0) AS ciendright,
        COALESCE(x_ci_in_start.n, 0) AS ci_in_start,
        COALESCE(x_ci_in_stop.n, 0) AS ci_in_stop,
        COALESCE(x_ci_start.n, 0) AS ci_start,
        COALESCE(x_ci_stop.n, 0) AS ci_stop,
        COALESCE(x_cytoband.n, 0) AS cytoband,
        COALESCE(x_cytoband_tax_id.n, 0) AS cytoband_tax_id,
        COALESCE(x_hgvs_name.n, 0) AS hgvs_name,
        COALESCE(x_plac_mutation_order.n, 0) AS plac_mutation_order,
        COALESCE(x_placement_rank.n, 0) AS placement_rank,
        COALESCE(x_placements_per_assembly.n, 0) AS placements_per_assembly,
        COALESCE(x_remap_best_within_cluster.n, 0) AS remap_best_within_cluster,
        COALESCE(x_remap_diff_chr.n, 0) AS remap_diff_chr,
        COALESCE(x_remap_failure_code.n, 0) AS remap_failure_code
    FROM sp
        LEFT JOIN (SELECT variant_placement_id, COUNT(alt_status) AS n FROM dgva.variant_placement_alt_status
                   GROUP BY variant_placement_id) x_alt_status ON x_alt_status.variant_placement_id = sp.variant_placement_id
        LEFT JOIN (SELECT variant_placement_id, COUNT(assembly_unit) AS n FROM dgva.variant_placement_assem_unit
                   GROUP BY variant_placement_id) x_assembly_unit ON x_assembly_unit.variant_placement_id = sp.variant_placement_id
        LEFT JOIN (SELECT variant_placement_id, COUNT(ciposleft) AS n FROM dgva.variant_placement_ciposleft
                   GROUP BY variant_placement_id) x_ciposleft ON x_ciposleft.variant_placement_id = sp.variant_placement_id
        LEFT JOIN (SELECT variant_placement_id, COUNT(ciposright) AS n FROM dgva.variant_placement_ciposright
                   GROUP BY variant_placement_id) x_ciposright ON x_ciposright.variant_placement_id = sp.variant_placement_id
        LEFT JOIN (SELECT variant_placement_id, COUNT(ciendleft) AS n FROM dgva.variant_placement_ciendleft
                   GROUP BY variant_placement_id) x_ciendleft ON x_ciendleft.variant_placement_id = sp.variant_placement_id
        LEFT JOIN (SELECT variant_placement_id, COUNT(ciendright) AS n FROM dgva.variant_placement_ciendright
                   GROUP BY variant_placement_id) x_ciendright ON x_ciendright.variant_placement_id = sp.variant_placement_id
        LEFT JOIN (SELECT variant_placement_id, COUNT(ci_plus) AS n FROM dgva.variant_placement_ci_in_start
                   GROUP BY variant_placement_id) x_ci_in_start ON x_ci_in_start.variant_placement_id = sp.variant_placement_id
        LEFT JOIN (SELECT variant_placement_id, COUNT(ci_plus) AS n FROM dgva.variant_placement_ci_in_stop
                   GROUP BY variant_placement_id) x_ci_in_stop ON x_ci_in_stop.variant_placement_id = sp.variant_placement_id
        LEFT JOIN (SELECT variant_placement_id, COUNT(ci_plus) AS n FROM dgva.variant_placement_ci_start
                   GROUP BY variant_placement_id) x_ci_start ON x_ci_start.variant_placement_id = sp.variant_placement_id
        LEFT JOIN (SELECT variant_placement_id, COUNT(ci_plus) AS n FROM dgva.variant_placement_ci_stop
                   GROUP BY variant_placement_id) x_ci_stop ON x_ci_stop.variant_placement_id = sp.variant_placement_id
        LEFT JOIN (SELECT variant_placement_id, COUNT(cytoband) AS n FROM dgva.variant_placement_cytoband
                   GROUP BY variant_placement_id) x_cytoband ON x_cytoband.variant_placement_id = sp.variant_placement_id
        LEFT JOIN (SELECT variant_placement_id, COUNT(cytoband_tax_id) AS n FROM dgva.variant_placement_cytoband
                   GROUP BY variant_placement_id) x_cytoband_tax_id ON x_cytoband_tax_id.variant_placement_id = sp.variant_placement_id
        LEFT JOIN (SELECT variant_placement_id, COUNT(hgvs_name) AS n FROM dgva.variant_placement_hgvs
                   GROUP BY variant_placement_id) x_hgvs_name ON x_hgvs_name.variant_placement_id = sp.variant_placement_id
        LEFT JOIN (SELECT variant_placement_id, COUNT(mutation_order) AS n FROM dgva.var_plac_mutation_order
                   GROUP BY variant_placement_id) x_plac_mutation_order ON x_plac_mutation_order.variant_placement_id = sp.variant_placement_id
        LEFT JOIN (SELECT variant_placement_id, COUNT(placement_rank) AS n FROM dgva.var_plac_placement_rank
                   GROUP BY variant_placement_id) x_placement_rank ON x_placement_rank.variant_placement_id = sp.variant_placement_id
        LEFT JOIN (SELECT variant_placement_id, COUNT(placements_per_assembly) AS n FROM dgva.var_plac_placs_per_assem
                   GROUP BY variant_placement_id) x_placements_per_assembly ON x_placements_per_assembly.variant_placement_id = sp.variant_placement_id
        LEFT JOIN (SELECT variant_placement_id, COUNT(remap_best_within_cluster) AS n FROM dgva.var_plac_remap_best_clust
                   GROUP BY variant_placement_id) x_remap_best_within_cluster ON x_remap_best_within_cluster.variant_placement_id = sp.variant_placement_id
        LEFT JOIN (SELECT variant_placement_id, COUNT(remap_diff_chr) AS n FROM dgva.var_plac_remap_diff_chr
                   GROUP BY variant_placement_id) x_remap_diff_chr ON x_remap_diff_chr.variant_placement_id = sp.variant_placement_id
        LEFT JOIN (SELECT variant_placement_id, COUNT(remap_failure_code) AS n FROM dgva.var_plac_remap_failure_code
                   GROUP BY variant_placement_id) x_remap_failure_code ON x_remap_failure_code.variant_placement_id = sp.variant_placement_id
),
per_variant AS (
    SELECT study_accession, region_or_call, variant_study_id, COUNT(*) AS placements,
        SUM(chr_accession) AS chr_accession_v, MAX(CASE WHEN chr_accession > 0 THEN 1 ELSE 0 END) AS chr_accession_var,
        SUM(alt_status) AS alt_status_v, MAX(CASE WHEN alt_status > 0 THEN 1 ELSE 0 END) AS alt_status_var,
        SUM(assembly_unit) AS assembly_unit_v, MAX(CASE WHEN assembly_unit > 0 THEN 1 ELSE 0 END) AS assembly_unit_var,
        SUM(ciposleft) AS ciposleft_v, MAX(CASE WHEN ciposleft > 0 THEN 1 ELSE 0 END) AS ciposleft_var,
        SUM(ciposright) AS ciposright_v, MAX(CASE WHEN ciposright > 0 THEN 1 ELSE 0 END) AS ciposright_var,
        SUM(ciendleft) AS ciendleft_v, MAX(CASE WHEN ciendleft > 0 THEN 1 ELSE 0 END) AS ciendleft_var,
        SUM(ciendright) AS ciendright_v, MAX(CASE WHEN ciendright > 0 THEN 1 ELSE 0 END) AS ciendright_var,
        SUM(ci_in_start) AS ci_in_start_v, MAX(CASE WHEN ci_in_start > 0 THEN 1 ELSE 0 END) AS ci_in_start_var,
        SUM(ci_in_stop) AS ci_in_stop_v, MAX(CASE WHEN ci_in_stop > 0 THEN 1 ELSE 0 END) AS ci_in_stop_var,
        SUM(ci_start) AS ci_start_v, MAX(CASE WHEN ci_start > 0 THEN 1 ELSE 0 END) AS ci_start_var,
        SUM(ci_stop) AS ci_stop_v, MAX(CASE WHEN ci_stop > 0 THEN 1 ELSE 0 END) AS ci_stop_var,
        SUM(cytoband) AS cytoband_v, MAX(CASE WHEN cytoband > 0 THEN 1 ELSE 0 END) AS cytoband_var,
        SUM(cytoband_tax_id) AS cytoband_tax_id_v, MAX(CASE WHEN cytoband_tax_id > 0 THEN 1 ELSE 0 END) AS cytoband_tax_id_var,
        SUM(hgvs_name) AS hgvs_name_v, MAX(CASE WHEN hgvs_name > 0 THEN 1 ELSE 0 END) AS hgvs_name_var,
        SUM(plac_mutation_order) AS plac_mutation_order_v, MAX(CASE WHEN plac_mutation_order > 0 THEN 1 ELSE 0 END) AS plac_mutation_order_var,
        SUM(placement_rank) AS placement_rank_v, MAX(CASE WHEN placement_rank > 0 THEN 1 ELSE 0 END) AS placement_rank_var,
        SUM(placements_per_assembly) AS placements_per_assembly_v, MAX(CASE WHEN placements_per_assembly > 0 THEN 1 ELSE 0 END) AS placements_per_assembly_var,
        SUM(remap_best_within_cluster) AS remap_best_within_cluster_v, MAX(CASE WHEN remap_best_within_cluster > 0 THEN 1 ELSE 0 END) AS remap_best_within_cluster_var,
        SUM(remap_diff_chr) AS remap_diff_chr_v, MAX(CASE WHEN remap_diff_chr > 0 THEN 1 ELSE 0 END) AS remap_diff_chr_var,
        SUM(remap_failure_code) AS remap_failure_code_v, MAX(CASE WHEN remap_failure_code > 0 THEN 1 ELSE 0 END) AS remap_failure_code_var
    FROM per_placement
    GROUP BY study_accession, region_or_call, variant_study_id
)
SELECT study_accession, region_or_call, COUNT(*) AS variants, SUM(placements) AS placements,
    SUM(chr_accession_v) AS chr_accession_v, SUM(chr_accession_var) AS chr_accession_var,
    SUM(alt_status_v) AS alt_status_v, SUM(alt_status_var) AS alt_status_var,
    SUM(assembly_unit_v) AS assembly_unit_v, SUM(assembly_unit_var) AS assembly_unit_var,
    SUM(ciposleft_v) AS ciposleft_v, SUM(ciposleft_var) AS ciposleft_var,
    SUM(ciposright_v) AS ciposright_v, SUM(ciposright_var) AS ciposright_var,
    SUM(ciendleft_v) AS ciendleft_v, SUM(ciendleft_var) AS ciendleft_var,
    SUM(ciendright_v) AS ciendright_v, SUM(ciendright_var) AS ciendright_var,
    SUM(ci_in_start_v) AS ci_in_start_v, SUM(ci_in_start_var) AS ci_in_start_var,
    SUM(ci_in_stop_v) AS ci_in_stop_v, SUM(ci_in_stop_var) AS ci_in_stop_var,
    SUM(ci_start_v) AS ci_start_v, SUM(ci_start_var) AS ci_start_var,
    SUM(ci_stop_v) AS ci_stop_v, SUM(ci_stop_var) AS ci_stop_var,
    SUM(cytoband_v) AS cytoband_v, SUM(cytoband_var) AS cytoband_var,
    SUM(cytoband_tax_id_v) AS cytoband_tax_id_v, SUM(cytoband_tax_id_var) AS cytoband_tax_id_var,
    SUM(hgvs_name_v) AS hgvs_name_v, SUM(hgvs_name_var) AS hgvs_name_var,
    SUM(plac_mutation_order_v) AS plac_mutation_order_v, SUM(plac_mutation_order_var) AS plac_mutation_order_var,
    SUM(placement_rank_v) AS placement_rank_v, SUM(placement_rank_var) AS placement_rank_var,
    SUM(placements_per_assembly_v) AS placements_per_assembly_v, SUM(placements_per_assembly_var) AS placements_per_assembly_var,
    SUM(remap_best_within_cluster_v) AS remap_best_within_cluster_v, SUM(remap_best_within_cluster_var) AS remap_best_within_cluster_var,
    SUM(remap_diff_chr_v) AS remap_diff_chr_v, SUM(remap_diff_chr_var) AS remap_diff_chr_var,
    SUM(remap_failure_code_v) AS remap_failure_code_v, SUM(remap_failure_code_var) AS remap_failure_code_var
FROM per_variant
GROUP BY study_accession, region_or_call
ORDER BY study_accession, region_or_call
