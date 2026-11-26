-- One row per study: whether it is listed in REMAPED_STUDIES and in GRCH38_STUDIES, the number of remapped
-- Exchange files stored in STUDY_FILE, and the types of all its stored files.
-- Binds: none.
SELECT
    s.study_accession,
    s.display_name,
    CASE WHEN EXISTS (SELECT 1 FROM dgva.remaped_studies x WHERE x.study_accession = s.study_accession)
         THEN 'Y' ELSE 'N' END AS in_remaped_studies,
    CASE WHEN EXISTS (SELECT 1 FROM dgva.grch38_studies x WHERE x.study_accession = s.study_accession)
         THEN 'Y' ELSE 'N' END AS in_grch38_studies,
    (SELECT COUNT(*) FROM dgva.study_file f WHERE f.study_accession = s.study_accession
        AND f.file_type LIKE 'Rem%') AS remapped_exchange_files,
    (SELECT LISTAGG(f.file_type, ',') WITHIN GROUP (ORDER BY f.file_type) FROM dgva.study_file f
        WHERE f.study_accession = s.study_accession) AS file_types
FROM dgva.dgva_study s
ORDER BY s.study_accession
