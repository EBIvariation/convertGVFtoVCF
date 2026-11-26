-- The files stored in STUDY_FILE for one study, with their type, date, checksum, size, number of 2000-byte chunks
-- and first bytes. FIRST_BYTES shows how the content is stored: b'<?xm' is plain XML, b'\x1f\x8b...' is gzip.
-- Binds: study (one study accession).
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
