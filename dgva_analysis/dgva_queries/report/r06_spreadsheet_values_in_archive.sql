-- Counts of values in the archive for: ALLELE_FREQUENCY rows of estd176; non-null REFERENCE_COPY_NUMBER in
-- estd206, estd221, estd222, estd224, estd228, estd232 and estd233; SUBJECT_FAMILY_RELATIONSHIP rows of estd195;
-- non-null SUBJECT_AGE in estd210, estd225 and estd226; subjects of estd222, and those with a SUBJECT_SEX.
-- Binds: none.
SELECT 'estd176 ALLELE_FREQUENCY rows' AS what, COUNT(*) AS n
FROM dgva.variant_study vs JOIN dgva.allele_frequency a ON a.variant_call_accession = vs.variant_accession
WHERE vs.study_accession = 'estd176'
UNION ALL
SELECT 'REFERENCE_COPY_NUMBER filled, ' || vs.study_accession, COUNT(c.reference_copy_number)
FROM dgva.variant_study vs JOIN dgva.variant_call_copy_number c ON c.variant_call_accession = vs.variant_accession
WHERE vs.study_accession IN ('estd206', 'estd221', 'estd222', 'estd224', 'estd228', 'estd232', 'estd233')
GROUP BY vs.study_accession
UNION ALL
SELECT 'estd195 SUBJECT_FAMILY_RELATIONSHIP rows', COUNT(*)
FROM dgva.dgva_subject s JOIN dgva.subject_family_relationship f ON f.subject_id = s.subject_id
WHERE s.study_accession = 'estd195'
UNION ALL
SELECT 'SUBJECT_AGE filled, ' || s.study_accession, COUNT(s.subject_age)
FROM dgva.dgva_subject s WHERE s.study_accession IN ('estd210', 'estd225', 'estd226') GROUP BY s.study_accession
UNION ALL
SELECT 'estd222 subjects', COUNT(*) FROM dgva.dgva_subject s WHERE s.study_accession = 'estd222'
UNION ALL
SELECT 'estd222 subjects with SUBJECT_SEX', COUNT(s.subject_sex) FROM dgva.dgva_subject s WHERE s.study_accession = 'estd222'
