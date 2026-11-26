-- Sessions opened by ReadOnlyOracleClient (MODULE 'convertGVFtoVCF-readonly'), with the statement each is running.
-- Needs SELECT on V$SESSION.
-- Binds: none.
SELECT s.sid, s.serial#, s.status, s.module, s.logon_time, s.last_call_et, s.sql_id, s.event
FROM v$session s
WHERE s.module = 'convertGVFtoVCF-readonly'
ORDER BY s.logon_time
