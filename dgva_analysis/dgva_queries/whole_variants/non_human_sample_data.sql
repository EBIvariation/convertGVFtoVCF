-- Studies whose taxonomy is not 9606, with their calls and how many have a sample, a sampleset, a zygosity and a
-- copy number.
-- Binds: none. Reads VARIANT_STUDY and the per-call sample tables.
WITH study_taxonomy AS (
    SELECT vs.study_accession, MIN(ga.assembly_taxonomy_id) AS taxonomy_id, MAX(ga.assembly_taxonomy_id) AS max_taxonomy_id
    FROM dgva.variant_study vs
        JOIN dgva.variant_placement p ON p.variant_study_id = vs.variant_study_id
        JOIN dgva.genome_assembly ga ON ga.assembly_accession = p.assembly_accession
    GROUP BY vs.study_accession
)
SELECT
    st.study_accession,
    st.taxonomy_id,
    st.max_taxonomy_id,
    COUNT(*) AS calls,
    SUM(CASE WHEN EXISTS (SELECT 1 FROM dgva.variant_call_sample s WHERE s.variant_call_accession = vs.variant_accession)
             THEN 1 ELSE 0 END) AS calls_with_sample,
    SUM(CASE WHEN EXISTS (SELECT 1 FROM dgva.variant_call_sampleset s WHERE s.variant_call_accession = vs.variant_accession)
             THEN 1 ELSE 0 END) AS calls_with_sampleset,
    SUM(CASE WHEN EXISTS (SELECT 1 FROM dgva.variant_call_zygosity z WHERE z.variant_call_accession = vs.variant_accession)
             THEN 1 ELSE 0 END) AS calls_with_zygosity,
    SUM(CASE WHEN EXISTS (SELECT 1 FROM dgva.variant_call_copy_number c WHERE c.variant_call_accession = vs.variant_accession)
             THEN 1 ELSE 0 END) AS calls_with_copy_number
FROM study_taxonomy st
    JOIN dgva.variant_study vs ON vs.study_accession = st.study_accession AND vs.region_or_call = 'call'
WHERE st.max_taxonomy_id <> 9606
GROUP BY st.study_accession, st.taxonomy_id, st.max_taxonomy_id
ORDER BY calls_with_sample DESC
