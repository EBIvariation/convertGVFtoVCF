-- Per study, for estd59, estd199, estd214 and estd219: calls, calls with a sample, and calls with a zygosity.
-- Binds: none.
SELECT vs.study_accession,
       COUNT(*) AS calls,
       COUNT(DISTINCT s.variant_call_accession) AS calls_with_sample,
       COUNT(DISTINCT z.variant_call_accession) AS calls_with_zygosity
FROM dgva.variant_study vs
    LEFT JOIN dgva.variant_call_sample s ON s.variant_call_accession = vs.variant_accession
    LEFT JOIN dgva.variant_call_zygosity z ON z.variant_call_accession = vs.variant_accession
WHERE vs.region_or_call = 'call'
    AND vs.study_accession IN ('estd59', 'estd199', 'estd214', 'estd219')
GROUP BY vs.study_accession
ORDER BY vs.study_accession
