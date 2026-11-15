-- Ticket 8.2. The estd214 and estd219 load files with copy numbers (ESTD214_CALLS_CN, ESTD219_CALLS_CN, and the AF
-- column of ESTD219_CALLS_CN): did the values reach VARIANT_CALL_COPY_NUMBER and ALLELE_FREQUENCY?
-- VARIANT_CALL_COPY_NUMBER.COPY_NUMBER is text; it is compared as a number when it looks like one.
-- Binds: none. Reads 15.8 M load rows.
WITH load_rows AS (
    SELECT 'ESTD214_CALLS_CN' AS load_table, accession, region_or_call, copy_number, CAST(NULL AS NUMBER) AS af
    FROM dgva.estd214_calls_cn
    UNION ALL
    SELECT 'ESTD219_CALLS_CN', accession, region_or_call, copy_number, af
    FROM dgva.estd219_calls_cn
)
SELECT
    l.load_table,
    l.region_or_call,
    COUNT(*) AS load_rows,
    SUM(CASE WHEN vs.variant_accession IS NOT NULL THEN 1 ELSE 0 END) AS accession_in_variant_study,
    SUM(CASE WHEN w.variant_accession IS NOT NULL THEN 1 ELSE 0 END) AS withdrawn,
    SUM(CASE WHEN l.copy_number IS NOT NULL THEN 1 ELSE 0 END) AS load_with_copy_number,
    SUM(CASE WHEN l.copy_number IS NOT NULL AND cn.variant_call_accession IS NOT NULL THEN 1 ELSE 0 END)
        AS copy_number_in_archive,
    SUM(CASE WHEN l.copy_number IS NOT NULL
              AND REGEXP_LIKE(TRIM(cn.copy_number), '^-?[0-9]+(\.[0-9]+)?$')
              AND TO_NUMBER(TRIM(cn.copy_number)) = l.copy_number THEN 1 ELSE 0 END) AS copy_number_equal,
    SUM(CASE WHEN l.af IS NOT NULL THEN 1 ELSE 0 END) AS load_with_af,
    SUM(CASE WHEN l.af IS NOT NULL AND af.variant_call_accession IS NOT NULL THEN 1 ELSE 0 END) AS af_in_archive,
    SUM(CASE WHEN l.af IS NOT NULL AND ABS(af.allele_frequency - l.af) < 1e-9 THEN 1 ELSE 0 END) AS af_equal
FROM load_rows l
    LEFT JOIN dgva.variant_study vs ON vs.variant_accession = l.accession
    LEFT JOIN dgva.variant_withdrawn w ON w.variant_accession = l.accession
    LEFT JOIN dgva.variant_call_copy_number cn ON cn.variant_call_accession = l.accession
    LEFT JOIN dgva.allele_frequency af ON af.variant_call_accession = l.accession
GROUP BY l.load_table, l.region_or_call
ORDER BY l.load_table, l.region_or_call
