-- Rows inserted into each archive table by each staging-to-archive transfer (STAGE_TO_ARCHIVE_ROWCOUNT), with the
-- study and the transfer and load IDs.
-- Binds: none. Small tables only.
SELECT sal.study_accession,
       sal.stage_to_archive_log_id,
       sal.stage_log_id,
       r.table_name,
       r.inserted_rowcount
FROM dgva.stage_to_archive_rowcount r
JOIN dgva.stage_to_archive_log sal ON sal.stage_to_archive_log_id = r.stage_to_archive_log_id
ORDER BY sal.study_accession, sal.stage_to_archive_log_id, r.table_name
