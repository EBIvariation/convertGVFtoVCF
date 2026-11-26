-- Variants of one study that are in VARIANT_STUDY_NEW only: their accession is in no VARIANT_STUDY row and not in
-- VARIANT_WITHDRAWN. One row per variant, with its kind, submitter ID and creation date.
-- Binds: study (one study accession). Reads the study's rows of VARIANT_STUDY_NEW (14.6 M rows in all).
SELECT n.variant_accession,
       n.region_or_call,
       n.submitter_variant_id,
       n.created_date
FROM dgva.variant_study_new n
WHERE n.study_accession = :study
    AND NOT EXISTS (SELECT 1 FROM dgva.variant_study a WHERE a.variant_accession = n.variant_accession)
    AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = n.variant_accession)
ORDER BY n.variant_accession
