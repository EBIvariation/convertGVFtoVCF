process PARSE_CREDENTIALS {
    tag "Parsing credentials" // prints a line

    input:
    val config_file

    output:
    tuple val(webin_user), val(webin_pass), emit: credentials

    exec:
    def target_file = file(config_file.toAbsolutePath())
    // Read from the TEST.config
    def yaml_data   = new org.yaml.snakeyaml.Yaml().load(target_file.text)
    webin_user      = yaml_data?.webin?.username
    webin_pass      = yaml_data?.webin?.password

    if ( !webin_user || !webin_pass ) {
        error "ERROR: Failed to extract credentials from ${config_file}"
    }
}