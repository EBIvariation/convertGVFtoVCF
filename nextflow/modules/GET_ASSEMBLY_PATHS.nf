process GET_ASSEMBLY_PATHS {
    tag "Getting assembly paths"

    input:
    path gvf_file
    
    output:
    tuple path(gvf_file), path("assembly.txt"), path("fasta.txt"), path("report.txt"), path("accession.txt")

    script:
    """
    export PYTHONPATH="${params.executable.convert_gvf.script_path}"

    if [[ "${gvf_file.name}" == *".p13."* ]]; then
        ASSEMBLY=\$(echo "${gvf_file.name}" | cut -d'.' -f3,4)
    else
        ASSEMBLY=\$(echo "${gvf_file.name}" | cut -d'.' -f3)
    fi

    CONVERT_CONFIG_PATH=\$(${params.executable.convert_gvf.interpreter} -c "import convert_gvf_to_vcf, os; print(os.path.join(os.path.dirname(convert_gvf_to_vcf.__file__), 'etc', 'config.yaml'))")

    ASSEMBLY_FASTA=\$(${params.executable.convert_gvf.interpreter} -c "import yaml, os; print(os.path.expandvars(yaml.safe_load(open('\$CONVERT_CONFIG_PATH'))['assembly_paths']['\$ASSEMBLY']))")
    ASSEMBLY_REPORT=\$(${params.executable.convert_gvf.interpreter} -c "import yaml, os; print(os.path.expandvars(yaml.safe_load(open('\$CONVERT_CONFIG_PATH'))['assembly_report_paths']['\$ASSEMBLY']))")
    ASSEMBLY_ACCESSION=\$(${params.executable.convert_gvf.interpreter} -c "import yaml, os; print(os.path.expandvars(yaml.safe_load(open('\$CONVERT_CONFIG_PATH'))['assembly_accession']['\$ASSEMBLY']))")
    
    echo "\$ASSEMBLY" > assembly.txt
    echo "\$ASSEMBLY_FASTA" > fasta.txt
    echo "\$ASSEMBLY_REPORT" > report.txt
    echo "\$ASSEMBLY_ACCESSION" > accession.txt
    """
}