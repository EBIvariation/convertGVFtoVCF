-- Pairs of VARIANT_REGION.EXCHANGE_ASSERTION_METHOD_ID and VAR_REG_ASSERT_METHOD.ASSERTION_METHOD, per study, with
-- the number of non-withdrawn regions for each pair.
-- Binds: none.
SELECT vs.study_accession, vr.exchange_assertion_method_id, am.assertion_method, COUNT(*) AS regions
FROM dgva.variant_region vr
    JOIN dgva.variant_study vs ON vs.variant_accession = vr.variant_region_accession
    JOIN dgva.var_reg_assert_method am ON am.var_reg_assert_method_id = vr.var_reg_assert_method_id
WHERE NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
GROUP BY vs.study_accession, vr.exchange_assertion_method_id, am.assertion_method
