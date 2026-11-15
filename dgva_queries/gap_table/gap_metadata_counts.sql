-- Ticket 7. Metadata candidate gap columns (reference: candidate gap list, rows 1-85, and row 112), counted per study:
-- rows with a non-null value, reached from the study through its experiments, samplesets, samples or subjects.
-- For a link table, the column is its second key: the count is the number of links.
-- STUDY_FILE.FILE_CONTENT is a BLOB: counted with IS NOT NULL like every other column.
-- Binds: none. Small tables only.
SELECT * FROM (
    SELECT 1 AS gap_row, 'DGVA_ANALYSIS' AS table_name, 'ANALYSIS_TYPE' AS column_name, e.study_accession AS study_accession,
        SUM(CASE WHEN t.analysis_type IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_experiment e JOIN dgva.experiment_analysis x ON x.experiment_id = e.experiment_id JOIN dgva.dgva_analysis t ON t.analysis_id = x.analysis_id
    GROUP BY e.study_accession
    UNION ALL
    SELECT 2 AS gap_row, 'DGVA_ANALYSIS' AS table_name, 'REFERENCE_VALUE' AS column_name, e.study_accession AS study_accession,
        SUM(CASE WHEN t.reference_value IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_experiment e JOIN dgva.experiment_analysis x ON x.experiment_id = e.experiment_id JOIN dgva.dgva_analysis t ON t.analysis_id = x.analysis_id
    GROUP BY e.study_accession
    UNION ALL
    SELECT 3 AS gap_row, 'DGVA_ANALYSIS' AS table_name, 'REFERENCE_TYPE' AS column_name, e.study_accession AS study_accession,
        SUM(CASE WHEN t.reference_type IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_experiment e JOIN dgva.experiment_analysis x ON x.experiment_id = e.experiment_id JOIN dgva.dgva_analysis t ON t.analysis_id = x.analysis_id
    GROUP BY e.study_accession
    UNION ALL
    SELECT 4 AS gap_row, 'DGVA_ANALYSIS' AS table_name, 'ANALYSIS_DESCRIPTION' AS column_name, e.study_accession AS study_accession,
        SUM(CASE WHEN t.analysis_description IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_experiment e JOIN dgva.experiment_analysis x ON x.experiment_id = e.experiment_id JOIN dgva.dgva_analysis t ON t.analysis_id = x.analysis_id
    GROUP BY e.study_accession
    UNION ALL
    SELECT 5 AS gap_row, 'DGVA_DETECTION' AS table_name, 'DETECTION_METHOD' AS column_name, e.study_accession AS study_accession,
        SUM(CASE WHEN t.detection_method IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_experiment e JOIN dgva.experiment_detection x ON x.experiment_id = e.experiment_id JOIN dgva.dgva_detection t ON t.detection_id = x.detection_id
    GROUP BY e.study_accession
    UNION ALL
    SELECT 6 AS gap_row, 'DGVA_DETECTION' AS table_name, 'DETECTION_DESCRIPTION' AS column_name, e.study_accession AS study_accession,
        SUM(CASE WHEN t.detection_description IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_experiment e JOIN dgva.experiment_detection x ON x.experiment_id = e.experiment_id JOIN dgva.dgva_detection t ON t.detection_id = x.detection_id
    GROUP BY e.study_accession
    UNION ALL
    SELECT 7 AS gap_row, 'DGVA_EXPERIMENT' AS table_name, 'SUBMITTER_EXPERIMENT_ID' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.submitter_experiment_id IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_experiment t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 8 AS gap_row, 'DGVA_EXPERIMENT' AS table_name, 'EXPERIMENT_TYPE' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.experiment_type IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_experiment t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 9 AS gap_row, 'DGVA_EXPERIMENT' AS table_name, 'EXPERIMENT_SITE' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.experiment_site IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_experiment t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 10 AS gap_row, 'DGVA_EXPERIMENT' AS table_name, 'EXPERIMENT_RESOLUTION' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.experiment_resolution IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_experiment t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 11 AS gap_row, 'DGVA_METHOD' AS table_name, 'METHOD_TYPE' AS column_name, e.study_accession AS study_accession,
        SUM(CASE WHEN t.method_type IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_experiment e JOIN dgva.experiment_method x ON x.experiment_id = e.experiment_id JOIN dgva.dgva_method t ON t.method_id = x.method_id
    GROUP BY e.study_accession
    UNION ALL
    SELECT 12 AS gap_row, 'DGVA_METHOD' AS table_name, 'METHOD_RESOLUTION' AS column_name, e.study_accession AS study_accession,
        SUM(CASE WHEN t.method_resolution IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_experiment e JOIN dgva.experiment_method x ON x.experiment_id = e.experiment_id JOIN dgva.dgva_method t ON t.method_id = x.method_id
    GROUP BY e.study_accession
    UNION ALL
    SELECT 13 AS gap_row, 'DGVA_METHOD' AS table_name, 'METHOD_DESCRIPTION' AS column_name, e.study_accession AS study_accession,
        SUM(CASE WHEN t.method_description IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_experiment e JOIN dgva.experiment_method x ON x.experiment_id = e.experiment_id JOIN dgva.dgva_method t ON t.method_id = x.method_id
    GROUP BY e.study_accession
    UNION ALL
    SELECT 14 AS gap_row, 'DGVA_PLATFORM' AS table_name, 'PLATFORM_NAME' AS column_name, e.study_accession AS study_accession,
        SUM(CASE WHEN t.platform_name IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_experiment e JOIN dgva.experiment_platform x ON x.experiment_id = e.experiment_id JOIN dgva.dgva_platform t ON t.platform_id = x.platform_id
    GROUP BY e.study_accession
    UNION ALL
    SELECT 15 AS gap_row, 'DGVA_SAMPLESET' AS table_name, 'SUBMITTER_SAMPLESET_ID' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.submitter_sampleset_id IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_sampleset t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 16 AS gap_row, 'DGVA_SAMPLESET' AS table_name, 'SAMPLESET_NAME' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.sampleset_name IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_sampleset t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 17 AS gap_row, 'DGVA_SAMPLESET' AS table_name, 'SAMPLESET_SIZE' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.sampleset_size IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_sampleset t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 18 AS gap_row, 'DGVA_SAMPLESET' AS table_name, 'SAMPLESET_TYPE' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.sampleset_type IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_sampleset t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 19 AS gap_row, 'DGVA_SAMPLESET' AS table_name, 'SAMPLESET_SEX' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.sampleset_sex IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_sampleset t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 20 AS gap_row, 'DGVA_SAMPLESET' AS table_name, 'SAMPLESET_POPULATION' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.sampleset_population IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_sampleset t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 21 AS gap_row, 'DGVA_SAMPLESET' AS table_name, 'SAMPLESET_DESCRIPTION' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.sampleset_description IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_sampleset t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 22 AS gap_row, 'DGVA_STUDY' AS table_name, 'BIOPROJECT_ACCESSION' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.bioproject_accession IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_study t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 23 AS gap_row, 'DGVA_STUDY' AS table_name, 'DBGAP_ID' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.dbgap_id IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_study t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 24 AS gap_row, 'DGVA_STUDY' AS table_name, 'EGA_ID' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.ega_id IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_study t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 25 AS gap_row, 'DGVA_STUDY' AS table_name, 'STUDY_URL' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.study_url IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_study t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 26 AS gap_row, 'DGVA_STUDY' AS table_name, 'HOLD_DATE' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.hold_date IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_study t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 27 AS gap_row, 'DGVA_STUDY' AS table_name, 'STUDY_COMPLETE' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.study_complete IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_study t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 28 AS gap_row, 'DGVA_STUDY' AS table_name, 'CREATION_DATE' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.creation_date IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_study t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 29 AS gap_row, 'DGVA_STUDY' AS table_name, 'SUBMITTER_VERSION' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.submitter_version IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_study t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 30 AS gap_row, 'DGVA_STUDY' AS table_name, 'IS_CURATED' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.is_curated IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_study t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 31 AS gap_row, 'DGVA_SUBMISSION' AS table_name, 'IS_TEST_DATA' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.is_test_data IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_submission t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 32 AS gap_row, 'DGVA_SUBMISSION' AS table_name, 'VLOAD_ID' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.vload_id IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_submission t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 33 AS gap_row, 'DGVA_SUBMISSION' AS table_name, 'EXCHANGE_SCHEMA_VERSION' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.exchange_schema_version IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_submission t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 34 AS gap_row, 'DGVA_SUBMISSION' AS table_name, 'STUDY_ID' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.study_id IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_submission t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 35 AS gap_row, 'DGVA_SUBMISSION' AS table_name, 'IS_INDIRECT_SUBMISSION' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.is_indirect_submission IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_submission t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 36 AS gap_row, 'DGVA_SUBMISSION' AS table_name, 'VALIDATION_SCHEMA_VERSION' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.validation_schema_version IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_submission t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 37 AS gap_row, 'EXPERIMENT_CURATION' AS table_name, 'CURATOR_NAME' AS column_name, e.study_accession AS study_accession,
        SUM(CASE WHEN t.curator_name IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_experiment e JOIN dgva.experiment_curation t ON t.experiment_id = e.experiment_id
    GROUP BY e.study_accession
    UNION ALL
    SELECT 38 AS gap_row, 'EXPERIMENT_CURATION' AS table_name, 'CURATOR_EMAIL' AS column_name, e.study_accession AS study_accession,
        SUM(CASE WHEN t.curator_email IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_experiment e JOIN dgva.experiment_curation t ON t.experiment_id = e.experiment_id
    GROUP BY e.study_accession
    UNION ALL
    SELECT 39 AS gap_row, 'EXPERIMENT_CURATION' AS table_name, 'CURATED_SET_NAME' AS column_name, e.study_accession AS study_accession,
        SUM(CASE WHEN t.curated_set_name IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_experiment e JOIN dgva.experiment_curation t ON t.experiment_id = e.experiment_id
    GROUP BY e.study_accession
    UNION ALL
    SELECT 40 AS gap_row, 'EXPERIMENT_CURATION' AS table_name, 'CURATED_SET_LINK' AS column_name, e.study_accession AS study_accession,
        SUM(CASE WHEN t.curated_set_link IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_experiment e JOIN dgva.experiment_curation t ON t.experiment_id = e.experiment_id
    GROUP BY e.study_accession
    UNION ALL
    SELECT 41 AS gap_row, 'EXPERIMENT_EXT_LINK' AS table_name, 'EXTERNAL_LINK_ID' AS column_name, e.study_accession AS study_accession,
        SUM(CASE WHEN t.external_link_id IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_experiment e JOIN dgva.experiment_ext_link t ON t.experiment_id = e.experiment_id
    GROUP BY e.study_accession
    UNION ALL
    SELECT 42 AS gap_row, 'EXPERIMENT_LIBRARY' AS table_name, 'LIBRARY_ABBREV' AS column_name, e.study_accession AS study_accession,
        SUM(CASE WHEN t.library_abbrev IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_experiment e JOIN dgva.experiment_library t ON t.experiment_id = e.experiment_id
    GROUP BY e.study_accession
    UNION ALL
    SELECT 43 AS gap_row, 'EXPERIMENT_LIBRARY' AS table_name, 'EXTERNAL_DB_ID' AS column_name, e.study_accession AS study_accession,
        SUM(CASE WHEN t.external_db_id IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_experiment e JOIN dgva.experiment_library t ON t.experiment_id = e.experiment_id
    GROUP BY e.study_accession
    UNION ALL
    SELECT 44 AS gap_row, 'EXPERIMENT_LIBRARY' AS table_name, 'EXTERNAL_DB_ABBREV_NAME' AS column_name, e.study_accession AS study_accession,
        SUM(CASE WHEN t.external_db_abbrev_name IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_experiment e JOIN dgva.experiment_library t ON t.experiment_id = e.experiment_id
    GROUP BY e.study_accession
    UNION ALL
    SELECT 45 AS gap_row, 'EXPERIMENT_LIBRARY' AS table_name, 'INVALID' AS column_name, e.study_accession AS study_accession,
        SUM(CASE WHEN t.invalid IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_experiment e JOIN dgva.experiment_library t ON t.experiment_id = e.experiment_id
    GROUP BY e.study_accession
    UNION ALL
    SELECT 46 AS gap_row, 'EXPERIMENT_LINK_URL' AS table_name, 'URL' AS column_name, e.study_accession AS study_accession,
        SUM(CASE WHEN t.url IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_experiment e JOIN dgva.experiment_link_url t ON t.experiment_id = e.experiment_id
    GROUP BY e.study_accession
    UNION ALL
    SELECT 47 AS gap_row, 'EXPERIMENT_PLATFORM' AS table_name, 'PLATFORM_ID' AS column_name, e.study_accession AS study_accession,
        SUM(CASE WHEN t.platform_id IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_experiment e JOIN dgva.experiment_platform t ON t.experiment_id = e.experiment_id
    GROUP BY e.study_accession
    UNION ALL
    SELECT 50 AS gap_row, 'LINKED_EXPERIMENT' AS table_name, 'LINKED_EXPERIMENT_ID' AS column_name, e.study_accession AS study_accession,
        SUM(CASE WHEN t.linked_experiment_id IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_experiment e JOIN dgva.linked_experiment t ON t.experiment_id = e.experiment_id
    GROUP BY e.study_accession
    UNION ALL
    SELECT 53 AS gap_row, 'PLATFORM_EXT_LINK' AS table_name, 'EXTERNAL_LINK_ID' AS column_name, e.study_accession AS study_accession,
        SUM(CASE WHEN t.external_link_id IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_experiment e JOIN dgva.experiment_platform x ON x.experiment_id = e.experiment_id JOIN dgva.platform_ext_link t ON t.platform_id = x.platform_id
    GROUP BY e.study_accession
    UNION ALL
    SELECT 54 AS gap_row, 'PUBMED_PUBLICATION' AS table_name, 'FIRST_AUTHOR' AS column_name, x.study_accession AS study_accession,
        SUM(CASE WHEN t.first_author IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.study_pubmed_publication x JOIN dgva.pubmed_publication t ON t.pubmed_id = x.pubmed_id
    GROUP BY x.study_accession
    UNION ALL
    SELECT 55 AS gap_row, 'PUBMED_PUBLICATION' AS table_name, 'ARTICLE_ABSTRACT' AS column_name, x.study_accession AS study_accession,
        SUM(CASE WHEN t.article_abstract IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.study_pubmed_publication x JOIN dgva.pubmed_publication t ON t.pubmed_id = x.pubmed_id
    GROUP BY x.study_accession
    UNION ALL
    SELECT 56 AS gap_row, 'SAMPLESET_EXT_LINK' AS table_name, 'EXTERNAL_LINK_ID' AS column_name, s.study_accession AS study_accession,
        SUM(CASE WHEN t.external_link_id IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_sampleset s JOIN dgva.sampleset_ext_link t ON t.sampleset_id = s.sampleset_id
    GROUP BY s.study_accession
    UNION ALL
    SELECT 57 AS gap_row, 'SAMPLESET_ORGANISM' AS table_name, 'SAMPLESET_TAXONOMY_ID' AS column_name, s.study_accession AS study_accession,
        SUM(CASE WHEN t.sampleset_taxonomy_id IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_sampleset s JOIN dgva.sampleset_organism t ON t.sampleset_id = s.sampleset_id
    GROUP BY s.study_accession
    UNION ALL
    SELECT 58 AS gap_row, 'SAMPLESET_PHENOTYPE' AS table_name, 'PHENOTYPE_CV_ID' AS column_name, s.study_accession AS study_accession,
        SUM(CASE WHEN t.phenotype_cv_id IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_sampleset s JOIN dgva.sampleset_phenotype t ON t.sampleset_id = s.sampleset_id
    GROUP BY s.study_accession
    UNION ALL
    SELECT 59 AS gap_row, 'SAMPLESET_PHENO_DESCRIPTION' AS table_name, 'SAMPLESET_PHENO_DESCRIPTION' AS column_name, s.study_accession AS study_accession,
        SUM(CASE WHEN t.sampleset_pheno_description IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_sampleset s JOIN dgva.sampleset_pheno_description t ON t.sampleset_id = s.sampleset_id
    GROUP BY s.study_accession
    UNION ALL
    SELECT 60 AS gap_row, 'SAMPLE_EXT_LINK' AS table_name, 'EXTERNAL_LINK_ID' AS column_name, s.study_accession AS study_accession,
        SUM(CASE WHEN t.external_link_id IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_sample s JOIN dgva.sample_ext_link t ON t.sample_id = s.sample_id
    GROUP BY s.study_accession
    UNION ALL
    SELECT 61 AS gap_row, 'STUDY_ALIAS' AS table_name, 'STUDY_ALIAS_NAME' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.study_alias_name IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.study_alias t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 62 AS gap_row, 'STUDY_COMMENT' AS table_name, 'STUDY_COMMENT' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.study_comment IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.study_comment t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 63 AS gap_row, 'STUDY_COMMENT' AS table_name, 'COMMENT_USER_NAME' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.comment_user_name IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.study_comment t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 64 AS gap_row, 'STUDY_COMMENT' AS table_name, 'COMMENT_TIMESTAMP' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.comment_timestamp IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.study_comment t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 65 AS gap_row, 'STUDY_CONTACT' AS table_name, 'CONTACT_EMAIL' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.contact_email IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.study_contact t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 66 AS gap_row, 'STUDY_CONTACT' AS table_name, 'CONTACT_PHONE' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.contact_phone IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.study_contact t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 67 AS gap_row, 'STUDY_CONTACT' AS table_name, 'AFFILIATION_NAME' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.affiliation_name IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.study_contact t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 68 AS gap_row, 'STUDY_CONTACT' AS table_name, 'AFFILIATION_ADDRESS' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.affiliation_address IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.study_contact t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 69 AS gap_row, 'STUDY_CONTACT' AS table_name, 'AFFILIATION_URL' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.affiliation_url IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.study_contact t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 70 AS gap_row, 'STUDY_FILE' AS table_name, 'FILE_NAME' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.file_name IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.study_file t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 71 AS gap_row, 'STUDY_FILE' AS table_name, 'FILE_TYPE' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.file_type IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.study_file t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 72 AS gap_row, 'STUDY_FILE' AS table_name, 'DATE_ADDED' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.date_added IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.study_file t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 73 AS gap_row, 'STUDY_FILE' AS table_name, 'FILE_CHECKSUM' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.file_checksum IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.study_file t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 74 AS gap_row, 'STUDY_FILE' AS table_name, 'FILE_DESCRIPTION' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.file_description IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.study_file t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 75 AS gap_row, 'STUDY_FILE' AS table_name, 'FILE_CONTENT' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.file_content IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.study_file t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 76 AS gap_row, 'STUDY_UPDATE' AS table_name, 'SUBMISSION_VERSION' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.submission_version IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.study_update t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 77 AS gap_row, 'STUDY_UPDATE' AS table_name, 'UPDATE_COMMENT' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.update_comment IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.study_update t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 78 AS gap_row, 'STUDY_UPDATE' AS table_name, 'NEW_FEATURE' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.new_feature IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.study_update t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 79 AS gap_row, 'STUDY_UPDATE' AS table_name, 'CORRECTION' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.correction IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.study_update t
    GROUP BY t.study_accession
    UNION ALL
    SELECT 80 AS gap_row, 'SUBJECT_COLLECTION' AS table_name, 'SUBJECT_COLLECTION_NAME' AS column_name, s.study_accession AS study_accession,
        SUM(CASE WHEN t.subject_collection_name IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_subject s JOIN dgva.subject_collection t ON t.subject_id = s.subject_id
    GROUP BY s.study_accession
    UNION ALL
    SELECT 81 AS gap_row, 'SUBJECT_EXT_LINK' AS table_name, 'EXTERNAL_LINK_ID' AS column_name, s.study_accession AS study_accession,
        SUM(CASE WHEN t.external_link_id IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_subject s JOIN dgva.subject_ext_link t ON t.subject_id = s.subject_id
    GROUP BY s.study_accession
    UNION ALL
    SELECT 82 AS gap_row, 'SUBJECT_FAMILY_RELATIONSHIP' AS table_name, 'FAMILY_RELATIONSHIP' AS column_name, s.study_accession AS study_accession,
        SUM(CASE WHEN t.family_relationship IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_subject s JOIN dgva.subject_family_relationship t ON t.subject_id = s.subject_id
    GROUP BY s.study_accession
    UNION ALL
    SELECT 83 AS gap_row, 'SUBJECT_PHENOTYPE' AS table_name, 'PHENOTYPE_CV_ID' AS column_name, s.study_accession AS study_accession,
        SUM(CASE WHEN t.phenotype_cv_id IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_subject s JOIN dgva.subject_phenotype t ON t.subject_id = s.subject_id
    GROUP BY s.study_accession
    UNION ALL
    SELECT 84 AS gap_row, 'SUBJECT_PHENOTYPE_DESCRIPTION' AS table_name, 'PHENOTYPE_DESCRIPTION' AS column_name, s.study_accession AS study_accession,
        SUM(CASE WHEN t.phenotype_description IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_subject s JOIN dgva.subject_phenotype_description t ON t.subject_id = s.subject_id
    GROUP BY s.study_accession
    UNION ALL
    SELECT 85 AS gap_row, 'SUBJECT_SAMPLE' AS table_name, 'SAMPLE_ID' AS column_name, s.study_accession AS study_accession,
        SUM(CASE WHEN t.sample_id IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.dgva_subject s JOIN dgva.subject_sample t ON t.subject_id = s.subject_id
    GROUP BY s.study_accession
    UNION ALL
    SELECT 112 AS gap_row, 'VARIANT_WITHDRAWN' AS table_name, 'WITHDRAWAL_REASON' AS column_name, t.study_accession AS study_accession,
        SUM(CASE WHEN t.withdrawal_reason IS NOT NULL THEN 1 ELSE 0 END) AS rows_with_value
    FROM dgva.variant_withdrawn t
    GROUP BY t.study_accession
)
WHERE rows_with_value > 0
ORDER BY gap_row, study_accession
