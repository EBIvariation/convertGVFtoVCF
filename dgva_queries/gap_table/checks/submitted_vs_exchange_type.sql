-- Ticket 7, open check (gap rows 86, 108): how often the submitter's variant type differs from the Exchange type.
-- Non-withdrawn variants (about 1.8 k more than the submitted ones: those with no submitted placement).
-- One row per (variant type, Exchange type, submitter type): variants and studies. No binds. Reads VARIANT_CALL and
-- VARIANT_REGION once each.
SELECT 'call' AS region_or_call, vc.exchange_variant_call_type AS exchange_type,
    vc.submitter_variant_call_type AS submitter_type,
    CASE WHEN vc.submitter_variant_call_type IS NULL THEN 'submitter null'
         WHEN vc.submitter_variant_call_type = vc.exchange_variant_call_type THEN 'identical'
         WHEN LOWER(REPLACE(REPLACE(vc.submitter_variant_call_type, '_', ' '), '-', ' '))
              = LOWER(REPLACE(REPLACE(vc.exchange_variant_call_type, '_', ' '), '-', ' ')) THEN 'case or separator only'
         ELSE 'different' END AS comparison,
    COUNT(*) AS variants, COUNT(DISTINCT vs.study_accession) AS studies, MIN(vs.study_accession) AS example_study
FROM dgva.variant_call vc
    JOIN dgva.variant_study vs ON vs.variant_accession = vc.variant_call_accession
WHERE NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
GROUP BY vc.exchange_variant_call_type, vc.submitter_variant_call_type
UNION ALL
SELECT 'region', vr.exchange_variant_region_type, vr.submitter_variant_region_type,
    CASE WHEN vr.submitter_variant_region_type IS NULL THEN 'submitter null'
         WHEN vr.submitter_variant_region_type = vr.exchange_variant_region_type THEN 'identical'
         WHEN LOWER(REPLACE(REPLACE(vr.submitter_variant_region_type, '_', ' '), '-', ' '))
              = LOWER(REPLACE(REPLACE(vr.exchange_variant_region_type, '_', ' '), '-', ' ')) THEN 'case or separator only'
         ELSE 'different' END,
    COUNT(*), COUNT(DISTINCT vs.study_accession), MIN(vs.study_accession)
FROM dgva.variant_region vr
    JOIN dgva.variant_study vs ON vs.variant_accession = vr.variant_region_accession
WHERE NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
GROUP BY vr.exchange_variant_region_type, vr.submitter_variant_region_type
