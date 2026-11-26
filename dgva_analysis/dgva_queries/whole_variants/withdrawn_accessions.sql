-- The withdrawn variants of one study, with their kind.
-- Binds: study (one study accession).
SELECT
    w.variant_accession,
    vs.region_or_call,
    w.withdrawal_reason
FROM dgva.variant_withdrawn w
    LEFT JOIN dgva.variant_study vs ON vs.variant_accession = w.variant_accession
WHERE w.study_accession = :study OR vs.study_accession = :study
ORDER BY w.variant_accession
