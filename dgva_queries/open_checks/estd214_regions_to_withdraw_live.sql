-- Ticket 8.2. The ESTD214_REGIONS_TO_WITHDRAW accessions that are still in VARIANT_STUDY and not withdrawn, with
-- their placements and linked calls. (The other 3,492 are in neither VARIANT_STUDY nor VARIANT_WITHDRAWN.)
-- Binds: none.
SELECT
    t.variant_region_accession,
    t.region_type,
    t.chrom,
    t.start_pos AS list_start,
    t.stop_pos AS list_stop,
    vs.study_accession,
    vs.created_date,
    (SELECT COUNT(*) FROM dgva.variant_placement p WHERE p.variant_study_id = vs.variant_study_id) AS placements,
    (SELECT COUNT(*) FROM dgva.linked_variant lv WHERE lv.variant_region_accession = t.variant_region_accession)
        AS linked_calls
FROM dgva.estd214_regions_to_withdraw t
    JOIN dgva.variant_study vs ON vs.variant_accession = t.variant_region_accession
WHERE NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = t.variant_region_accession)
