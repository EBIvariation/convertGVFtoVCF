-- VLOAD ID and display name of every study, one row per study and submission (DGVA_SUBMISSION). Studies with no
-- DGVA_SUBMISSION row are kept, with a null VLOAD_ID.
-- Binds: none. Small tables only.
SELECT s.study_accession,
       s.display_name,
       s.study_type,
       sub.submission_id,
       sub.vload_id,
       sub.is_test_data,
       sub.is_indirect_submission,
       sub.exchange_schema_version
FROM dgva.dgva_study s
LEFT JOIN dgva.dgva_submission sub ON sub.study_accession = s.study_accession
ORDER BY CASE WHEN s.study_accession LIKE 'estd%' THEN 1 ELSE 2 END,
         TO_NUMBER(REGEXP_SUBSTR(s.study_accession, '[0-9]+')),
         sub.submission_id
