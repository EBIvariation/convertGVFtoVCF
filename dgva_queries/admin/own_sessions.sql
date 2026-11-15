-- Sessions opened by the read-only helper (MODULE set by ReadOnlyOracleClient), with what they are running.
-- Used to check that a stopped query is no longer running on the server. Needs SELECT on V$SESSION.
-- Binds: none.
SELECT s.sid, s.serial#, s.status, s.module, s.logon_time, s.last_call_et, s.sql_id, s.event
FROM v$session s
WHERE s.module = 'convertGVFtoVCF-readonly'
ORDER BY s.logon_time
