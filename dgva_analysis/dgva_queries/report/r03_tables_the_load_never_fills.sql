-- Row counts of VARIANT_CALL_MODE_OF_INHERIT, VARIANT_MATE and VARIANT_EVENT.
-- Binds: none. Small tables only.
SELECT 'VARIANT_CALL_MODE_OF_INHERIT' AS table_name, COUNT(*) AS row_count FROM dgva.variant_call_mode_of_inherit
UNION ALL SELECT 'VARIANT_MATE', COUNT(*) FROM dgva.variant_mate
UNION ALL SELECT 'VARIANT_EVENT', COUNT(*) FROM dgva.variant_event
