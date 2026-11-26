-- Inserts, updates and deletes recorded per table (ALL_TAB_MODIFICATIONS) since its statistics were last gathered.
-- Oracle writes this data to the dictionary periodically, so the latest changes may be missing.
-- Binds: owner (the schema, for example DGVA).
SELECT table_name, partition_name, inserts, updates, deletes, truncated, timestamp
FROM all_tab_modifications
WHERE table_owner = :owner
ORDER BY timestamp DESC, table_name
