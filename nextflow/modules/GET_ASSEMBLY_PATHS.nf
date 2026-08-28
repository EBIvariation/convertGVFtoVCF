process GET_ASSEMBLY_PATHS {
    tag "Getting assembly paths"

    input:
    path gvf_file

    output:
    tuple path(gvf_file), path("assembly.txt"), path("fasta.txt"), path("report.txt"), path("accession.txt")

    script:

    def p = gvf_file.name.tokenize('.')
    def ASSEMBLY = gvf_file.name.contains(".p13.") ? "${p[2]}.${p[3]}" : p[2]

    def CONVERT_CONFIG_PATH = "${params.executable.convert_gvf.config}"

    def yaml_data = new org.yaml.snakeyaml.Yaml().load(new File(CONVERT_CONFIG_PATH).text)

    def ASSEMBLY_FASTA     = yaml_data?.assembly_paths?."${ASSEMBLY}"
    def ASSEMBLY_REPORT    = yaml_data?.assembly_report_paths?."${ASSEMBLY}"
    def ASSEMBLY_ACCESSION = yaml_data?.assembly_accession?."${ASSEMBLY}"

    """
    export REF_PATH="${params.REF_PATH}"
    echo "${ASSEMBLY}" > assembly.txt
    echo "${ASSEMBLY_FASTA}" > fasta.txt
    echo "${ASSEMBLY_REPORT}" > report.txt
    echo "${ASSEMBLY_ACCESSION}" > accession.txt
    """
}
