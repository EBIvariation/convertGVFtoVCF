process SUBMIT_TO_EVA {
    tag "Submitting to EVA after validation SUCCESS: ${log_file.baseName}"
    
    input:
    path log_file
    tuple val(webin_user), val(webin_pass)
    val submit_study_dir
    val metadata_json

    when:
    params.RUN_SUBMIT

    script:
    """
    source ${params.executable.eva_sub_cli.script_path}
    eva-sub-cli.py --submission_dir "${submit_study_dir}" \
		           --metadata_json "${metadata_json}" \
		           --tasks submit \
		           --username "${webin_user}" \
		           --password "${webin_pass}"
    """
}