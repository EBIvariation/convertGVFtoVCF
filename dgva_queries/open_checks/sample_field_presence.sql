-- Ticket 8.3. For one study: how much data the database holds for the fields the B bugs drop from the GVF, so that
-- the FTP files can be judged (a bug only loses data where the database has some).
-- Counts are on non-withdrawn variants and on submitted placements (ticket 2 rule).
--   B2 submitted sequences; B3 confidence intervals (both sets; CIEND left <> right); B5 subject phenotypes, subject
--   and sample links; B6 samples outside any sampleset; B7 validations; B9 mates; B12 HGVS names.
-- Binds: study (one study accession). One row.
WITH v AS (
    SELECT vs.variant_study_id, vs.variant_accession, vs.region_or_call
    FROM dgva.variant_study vs
    WHERE vs.study_accession = :study
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
),
p AS (
    SELECT p.variant_placement_id, v.region_or_call
    FROM v
        JOIN dgva.variant_placement p ON p.variant_study_id = v.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.variant_placement_remap rm ON rm.variant_placement_id = p.variant_placement_id
    WHERE rm.variant_placement_id IS NULL
        AND (m.variant_placement_method IS NULL
             OR (m.variant_placement_method LIKE 'Submitted%'
                 AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                        'Submitted on non-INSDC assembly')))
)
SELECT
    :study AS study_accession,
    (SELECT COUNT(*) FROM v WHERE region_or_call = 'call') AS calls,
    (SELECT COUNT(*) FROM v JOIN dgva.variant_call_sequence s ON s.variant_call_accession = v.variant_accession
     WHERE s.variant_call_seq IS NOT NULL) AS b2_calls_with_sequence,
    (SELECT COUNT(*) FROM p) AS submitted_placements,
    (SELECT COUNT(*) FROM p JOIN dgva.variant_placement_ci_start c ON c.variant_placement_id = p.variant_placement_id
     WHERE p.region_or_call = 'call') AS b3_call_placements_with_ci_start,
    (SELECT COUNT(*) FROM p JOIN dgva.variant_placement_ciposleft c ON c.variant_placement_id = p.variant_placement_id
     WHERE p.region_or_call = 'call') AS b3_call_placements_with_cipos,
    (SELECT COUNT(*) FROM p JOIN dgva.variant_placement_ci_start c ON c.variant_placement_id = p.variant_placement_id
     WHERE p.region_or_call = 'region') AS b3_region_placements_with_ci_start,
    (SELECT COUNT(*) FROM p JOIN dgva.variant_placement_ciendleft l ON l.variant_placement_id = p.variant_placement_id
         JOIN dgva.variant_placement_ciendright r ON r.variant_placement_id = p.variant_placement_id
     WHERE p.region_or_call = 'region' AND l.ciendleft <> r.ciendright) AS b3_regions_ciend_left_ne_right,
    (SELECT COUNT(*) FROM p JOIN dgva.variant_placement_hgvs h ON h.variant_placement_id = p.variant_placement_id)
        AS b12_placements_with_hgvs,
    (SELECT COUNT(*) FROM v JOIN dgva.variant_call_validation c ON c.variant_call_accession = v.variant_accession)
        AS b7_validation_rows,
    (SELECT COUNT(*) FROM v JOIN dgva.variant_mate m ON m.variant_accession = v.variant_accession
     WHERE m.mate_accession <> v.variant_accession) AS b9_variants_with_other_mate,
    (SELECT COUNT(*) FROM dgva.dgva_subject s WHERE s.study_accession = :study) AS subjects,
    (SELECT COUNT(DISTINCT sp.subject_id) FROM dgva.dgva_subject s
         JOIN dgva.subject_phenotype sp ON sp.subject_id = s.subject_id WHERE s.study_accession = :study)
        AS b5_subjects_with_phenotype,
    (SELECT COUNT(DISTINCT sl.subject_id) FROM dgva.dgva_subject s
         JOIN dgva.subject_ext_link sl ON sl.subject_id = s.subject_id WHERE s.study_accession = :study)
        AS b5_subjects_with_link,
    (SELECT COUNT(*) FROM dgva.dgva_sample s WHERE s.study_accession = :study) AS samples,
    (SELECT COUNT(DISTINCT ss.sample_id) FROM dgva.dgva_sample s
         JOIN dgva.subject_sample ss ON ss.sample_id = s.sample_id WHERE s.study_accession = :study)
        AS b5_samples_linked_to_subject,
    (SELECT COUNT(DISTINCT sl.sample_id) FROM dgva.dgva_sample s
         JOIN dgva.sample_ext_link sl ON sl.sample_id = s.sample_id WHERE s.study_accession = :study)
        AS b5_samples_with_link,
    (SELECT COUNT(*) FROM dgva.dgva_sample s WHERE s.study_accession = :study
         AND NOT EXISTS (SELECT 1 FROM dgva.sample_sampleset x WHERE x.sample_id = s.sample_id))
        AS b6_samples_outside_sampleset
FROM dual
