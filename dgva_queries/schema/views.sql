-- Views in the schema.
-- Binds: owner (DGVA)
SELECT view_name
FROM all_views
WHERE owner = :owner
ORDER BY view_name
