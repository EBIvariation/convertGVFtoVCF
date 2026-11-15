-- Ticket 6, FTP check. The withdrawn variants of a study, with their type, to compare with the accessions in the
-- FTP GVFs (a file exported before a withdrawal still holds the variant).
-- Binds: study (one study accession).
SELECT
    w.variant_accession,
    vs.region_or_call,
    w.withdrawal_reason
FROM dgva.variant_withdrawn w
    LEFT JOIN dgva.variant_study vs ON vs.variant_accession = w.variant_accession
WHERE w.study_accession = :study OR vs.study_accession = :study
ORDER BY w.variant_accession
