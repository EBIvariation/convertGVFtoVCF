process VALIDATE_SUBMISSION {
    tag "Validating submission data: ${input_study_name}"
    
    publishDir { "${params.output_dir}/validation_reports/${input_study_name}" }, mode: 'copy'

    input:
    val input_study_name
    val study_accession
    val trigger_token 

    output:
    path "**/*.json", emit: validation_manifests, optional: true
    path ".command.log", emit: validation_log

    script:
    def submit_study_dir = "${params.output_dir}/submission/${input_study_name}"
    def metadata_json = "${submit_study_dir}/eva_submission_${study_accession}.json"
    """
    export PYTHONPATH="${params.executable.eva_sub_cli.script_path}"
    ${params.executable.eva_sub_cli.interpreter} ${params.executable.eva_sub_cli.script_path}/eva-sub-cli.py \\
        --submission_dir "${submit_study_dir}" \\
        --metadata_json "${metadata_json}" \\
        --tasks validate \\
        --validation_tasks ${params.VALIDATION_TASKS}
    """
}
