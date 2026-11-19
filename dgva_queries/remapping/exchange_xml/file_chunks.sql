-- Ticket 2, follow-up. The content of one STUDY_FILE row, as 2000-byte chunks (2000 bytes is the largest RAW that
-- DBMS_LOB.SUBSTR returns in SQL). The helper writes each chunk as a Python bytes literal; find_variants_in_exchange.py
-- joins them back into the file.
-- Binds: study, file_name (exactly as listed by study_files.sql), first_chunk, last_chunk (1-based, inclusive).
-- For the whole file: first_chunk=1, last_chunk=CHUNKS from study_files.out. Split very large files into several
-- runs (for example 50,000 chunks = 100 MB each) and pass all the outputs to the script.
WITH f AS (
    SELECT file_content, DBMS_LOB.GETLENGTH(file_content) AS size_bytes
    FROM dgva.study_file
    WHERE study_accession = :study AND file_name = :file_name
),
n AS (
    SELECT TO_NUMBER(:first_chunk) + LEVEL - 1 AS chunk
    FROM dual
    CONNECT BY LEVEL <= TO_NUMBER(:last_chunk) - TO_NUMBER(:first_chunk) + 1
)
SELECT
    n.chunk,
    DBMS_LOB.SUBSTR(f.file_content, 2000, (n.chunk - 1) * 2000 + 1) AS bytes
FROM f
    CROSS JOIN n
WHERE (n.chunk - 1) * 2000 < f.size_bytes
ORDER BY n.chunk
