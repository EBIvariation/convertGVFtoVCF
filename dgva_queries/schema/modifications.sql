-- Inserts, updates and deletes recorded per table since its statistics were last gathered.
-- Oracle writes this monitoring data to the dictionary periodically, so the latest changes may not show yet.
-- Binds: owner (DGVA)
SELECT table_name, partition_name, inserts, updates, deletes, truncated, timestamp
FROM all_tab_modifications
WHERE table_owner = :owner
ORDER BY timestamp DESC, table_name
