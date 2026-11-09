-- Tables in the schema and their optimiser row counts (approximate: check LAST_ANALYZED).
-- Binds: owner (DGVA)
SELECT table_name, num_rows, last_analyzed
FROM all_tables
WHERE owner = :owner
ORDER BY table_name
