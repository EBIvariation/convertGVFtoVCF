-- Ticket 6, W1 detail. LINKED_VARIANT rows by what they link: the parent is always a region
-- (VARIANT_REGION_ACCESSION); the linked accession can be a call, another region, or unknown.
-- Also says whether both ends are in the same study, and whether the linked variant is withdrawn.
-- Binds: study (the parent region's study accession, or ALL).
SELECT
    rs.study_accession AS region_study,
    CASE WHEN lc.variant_call_accession IS NOT NULL THEN 'region_to_call'
         WHEN lr.variant_region_accession IS NOT NULL THEN 'region_to_region'
         ELSE 'region_to_unknown' END AS link_kind,
    CASE WHEN ls.study_accession IS NULL THEN 'not_in_variant_study'
         WHEN ls.study_accession = rs.study_accession THEN 'same_study'
         ELSE 'other_study' END AS linked_study,
    CASE WHEN w.variant_accession IS NULL THEN 'N' ELSE 'Y' END AS linked_withdrawn,
    COUNT(*) AS links
FROM dgva.linked_variant lv
    LEFT JOIN dgva.variant_study rs ON rs.variant_accession = lv.variant_region_accession
    LEFT JOIN dgva.variant_call lc ON lc.variant_call_accession = lv.linked_variant_accession
    LEFT JOIN dgva.variant_region lr ON lr.variant_region_accession = lv.linked_variant_accession
    LEFT JOIN dgva.variant_study ls ON ls.variant_accession = lv.linked_variant_accession
    LEFT JOIN dgva.variant_withdrawn w ON w.variant_accession = lv.linked_variant_accession
WHERE :study = 'ALL' OR rs.study_accession = :study
GROUP BY rs.study_accession,
    CASE WHEN lc.variant_call_accession IS NOT NULL THEN 'region_to_call'
         WHEN lr.variant_region_accession IS NOT NULL THEN 'region_to_region'
         ELSE 'region_to_unknown' END,
    CASE WHEN ls.study_accession IS NULL THEN 'not_in_variant_study'
         WHEN ls.study_accession = rs.study_accession THEN 'same_study'
         ELSE 'other_study' END,
    CASE WHEN w.variant_accession IS NULL THEN 'N' ELSE 'Y' END
ORDER BY region_study, link_kind, linked_study, linked_withdrawn
