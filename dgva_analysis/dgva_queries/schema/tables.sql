-- Tables of a schema with their optimiser row counts (NUM_ROWS) and LAST_ANALYZED date.
-- Binds: owner (the schema, for example DGVA).
SELECT table_name, num_rows, last_analyzed
FROM all_tables
WHERE owner = :owner
ORDER BY table_name
