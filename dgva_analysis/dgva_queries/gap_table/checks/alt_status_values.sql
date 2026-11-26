-- Values of VARIANT_PLACEMENT_ALT_STATUS.ALT_STATUS, with each placement's study, assembly, chromosome or contig,
-- placement method and assembly unit.
-- Binds: none. Reads about 10 k rows.
SELECT a.alt_status, vs.study_accession, vs.region_or_call, p.assembly_accession,
    CASE WHEN p.chromosome_name IS NOT NULL THEN 'chromosome' WHEN c.contig_id IS NOT NULL THEN 'contig'
         ELSE 'none' END AS sequence_kind,
    NVL(m.variant_placement_method, '(null)') AS placement_method, NVL(u.assembly_unit, '(null)') AS assembly_unit,
    COUNT(*) AS placements, MIN(vs.variant_accession) AS example_variant,
    MIN(NVL(p.chromosome_name, c.contig_id)) AS example_sequence
FROM dgva.variant_placement_alt_status a
    JOIN dgva.variant_placement p ON p.variant_placement_id = a.variant_placement_id
    JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
    LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_contig c ON c.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_assem_unit u ON u.variant_placement_id = p.variant_placement_id
GROUP BY a.alt_status, vs.study_accession, vs.region_or_call, p.assembly_accession,
    CASE WHEN p.chromosome_name IS NOT NULL THEN 'chromosome' WHEN c.contig_id IS NOT NULL THEN 'contig'
         ELSE 'none' END,
    NVL(m.variant_placement_method, '(null)'), NVL(u.assembly_unit, '(null)')
