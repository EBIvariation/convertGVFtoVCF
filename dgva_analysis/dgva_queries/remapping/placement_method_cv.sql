-- The placement method vocabulary (VARIANT_PLACEMENT_METHOD_CV), with each method's METHOD_CLASS.
-- Binds: none.
SELECT variant_placement_method, method_class, is_current, added_exchange_version, removed_exchange_version
FROM dgva.variant_placement_method_cv
ORDER BY method_class, variant_placement_method
