-- Columns of every table, to compare with DGVA_DDL.sql.
-- Binds: owner (DGVA)
SELECT table_name, column_id, column_name, data_type, data_length, nullable
FROM all_tab_columns
WHERE owner = :owner
ORDER BY table_name, column_id
