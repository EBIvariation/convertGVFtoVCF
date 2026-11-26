-- Values of EXTERNAL_LINK.LINKED_OBJECT_TYPE, with the external database names and the link tables (calls,
-- samples, experiments, platforms, samplesets, subjects) that point at each link.
-- Binds: none. Reads the whole EXTERNAL_LINK table (about 6.9 M rows).
SELECT el.linked_object_type, el.external_db_name,
    CASE WHEN EXISTS (SELECT 1 FROM dgva.variant_call_ext_link x WHERE x.external_link_id = el.external_link_id) THEN 'call'
         WHEN EXISTS (SELECT 1 FROM dgva.sample_ext_link x WHERE x.external_link_id = el.external_link_id) THEN 'sample'
         WHEN EXISTS (SELECT 1 FROM dgva.experiment_ext_link x WHERE x.external_link_id = el.external_link_id) THEN 'experiment'
         WHEN EXISTS (SELECT 1 FROM dgva.platform_ext_link x WHERE x.external_link_id = el.external_link_id) THEN 'platform'
         WHEN EXISTS (SELECT 1 FROM dgva.sampleset_ext_link x WHERE x.external_link_id = el.external_link_id) THEN 'sampleset'
         WHEN EXISTS (SELECT 1 FROM dgva.subject_ext_link x WHERE x.external_link_id = el.external_link_id) THEN 'subject'
         ELSE 'none of these' END AS linked_from,
    COUNT(*) AS links, SUM(NVL(el.invalid, 0)) AS invalid_links, MIN(el.external_db_id) AS example_id
FROM dgva.external_link el
GROUP BY el.linked_object_type, el.external_db_name,
    CASE WHEN EXISTS (SELECT 1 FROM dgva.variant_call_ext_link x WHERE x.external_link_id = el.external_link_id) THEN 'call'
         WHEN EXISTS (SELECT 1 FROM dgva.sample_ext_link x WHERE x.external_link_id = el.external_link_id) THEN 'sample'
         WHEN EXISTS (SELECT 1 FROM dgva.experiment_ext_link x WHERE x.external_link_id = el.external_link_id) THEN 'experiment'
         WHEN EXISTS (SELECT 1 FROM dgva.platform_ext_link x WHERE x.external_link_id = el.external_link_id) THEN 'platform'
         WHEN EXISTS (SELECT 1 FROM dgva.sampleset_ext_link x WHERE x.external_link_id = el.external_link_id) THEN 'sampleset'
         WHEN EXISTS (SELECT 1 FROM dgva.subject_ext_link x WHERE x.external_link_id = el.external_link_id) THEN 'subject'
         ELSE 'none of these' END
