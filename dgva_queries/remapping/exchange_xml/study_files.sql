-- Ticket 2, follow-up. The files DGVa stored for a study in STUDY_FILE, with their size and first bytes.
-- FIRST_BYTES tells how the content is stored: b'<?xm' is plain XML, b'\x1f\x8b...' is gzip.
-- SIZE_BYTES gives the number of 2000-byte chunks to request with file_chunks.sql: CEIL(SIZE_BYTES / 2000).
-- Binds: study (a study accession, for example nstd35).
SELECT
    f.study_accession,
    f.file_name,
    f.file_type,
    f.date_added,
    f.file_checksum,
    DBMS_LOB.GETLENGTH(f.file_content) AS size_bytes,
    CEIL(DBMS_LOB.GETLENGTH(f.file_content) / 2000) AS chunks,
    DBMS_LOB.SUBSTR(f.file_content, 8, 1) AS first_bytes
FROM dgva.study_file f
WHERE f.study_accession = :study
ORDER BY f.file_type, f.date_added, f.file_name
