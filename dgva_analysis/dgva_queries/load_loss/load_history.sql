-- Every load recorded in STAGE_LOG, with the loaded file name, loading date, code and schema versions, submission
-- details, and the STAGE_TO_ARCHIVE_LOG transfers that followed. Loading user names are not selected.
-- Binds: none. Small tables only.
SELECT sl.study_accession,
       sl.stage_log_id,
       sl.loaded_file_name,
       sl.loading_date,
       sl.loading_code_version,
       sl.schema_version,
       sl.submission_file_name,
       sl.submission_ticket,
       sl.submission_version,
       sl.stage_log_comment,
       sal.stage_to_archive_log_id,
       sal.plsql_version,
       sal.stage_to_archive_log_comment
FROM dgva.stage_log sl
LEFT JOIN dgva.stage_to_archive_log sal ON sal.stage_log_id = sl.stage_log_id
ORDER BY sl.study_accession, sl.stage_log_id, sal.stage_to_archive_log_id
