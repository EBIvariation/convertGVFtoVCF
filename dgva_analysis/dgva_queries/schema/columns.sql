-- Columns of every table of a schema, in column order.
-- Binds: owner (the schema, for example DGVA).
SELECT table_name, column_id, column_name, data_type, data_length, nullable
FROM all_tab_columns
WHERE owner = :owner
ORDER BY table_name, column_id
