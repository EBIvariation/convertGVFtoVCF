-- Ticket 7, open check (gap row 107): is VARIANT_REGION.EXCHANGE_ASSERTION_METHOD_ID just the submission's own number
-- for the assertion method already written to the GVF (VAR_REG_ASSERT_METHOD.ASSERTION_METHOD)?
-- generate_exchange_archive.pl writes it back as <ASSERTION_METHOD assertion_method_id=...>.
-- One row per study and (Exchange ID, method) pair, with region counts, for non-withdrawn regions. No binds.
SELECT vs.study_accession, vr.exchange_assertion_method_id, am.assertion_method, COUNT(*) AS regions
FROM dgva.variant_region vr
    JOIN dgva.variant_study vs ON vs.variant_accession = vr.variant_region_accession
    JOIN dgva.var_reg_assert_method am ON am.var_reg_assert_method_id = vr.var_reg_assert_method_id
WHERE NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
GROUP BY vs.study_accession, vr.exchange_assertion_method_id, am.assertion_method
