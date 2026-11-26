-- Views of a schema.
-- Binds: owner (the schema, for example DGVA).
SELECT view_name
FROM all_views
WHERE owner = :owner
ORDER BY view_name
