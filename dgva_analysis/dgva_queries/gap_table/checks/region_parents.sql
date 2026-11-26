-- Region-to-region links in LINKED_VARIANT, one row per link: study, child and parent region, whether the child
-- has a submitter ID, whether the child or the parent is withdrawn, and the number of parent regions of the child.
-- Binds: none.
SELECT vs.study_accession, lv.linked_variant_accession AS child_region, lv.variant_region_accession AS parent_region,
    CASE WHEN vs.submitter_variant_id IS NULL THEN 'N' ELSE 'Y' END AS child_has_submitter_id,
    CASE WHEN EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = lv.linked_variant_accession)
         THEN 'Y' ELSE 'N' END AS child_withdrawn,
    CASE WHEN EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = lv.variant_region_accession)
         THEN 'Y' ELSE 'N' END AS parent_withdrawn,
    (SELECT COUNT(*) FROM dgva.linked_variant l2 WHERE l2.linked_variant_accession = lv.linked_variant_accession
        AND l2.variant_region_accession IN (SELECT variant_region_accession FROM dgva.variant_region)) AS parents_of_child
FROM dgva.linked_variant lv
    JOIN dgva.variant_region child ON child.variant_region_accession = lv.linked_variant_accession
    JOIN dgva.variant_study vs ON vs.variant_accession = lv.linked_variant_accession
