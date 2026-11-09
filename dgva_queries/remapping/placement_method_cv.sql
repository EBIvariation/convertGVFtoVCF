-- Ticket 2. The placement method vocabulary, including DGVa's own METHOD_CLASS for each method.
-- Binds: none
SELECT variant_placement_method, method_class, is_current, added_exchange_version, removed_exchange_version
FROM dgva.variant_placement_method_cv
ORDER BY method_class, variant_placement_method
