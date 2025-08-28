process SVANNA_PRIORITIZE {
    tag "$meta.id"
    label 'process_high'
    
    container "docker.io/nourmahfel1/svanna:latest"

    input:
    tuple val(meta), path(vcf)
    path(data_directory)
    val(hpo_terms)

    output:
    tuple val(meta), path("*.html"), emit: html_report, optional: true
    tuple val(meta), path("*.csv"), emit: csv_results, optional: true  
    tuple val(meta), path("*.vcf.gz"), emit: vcf_results, optional: true
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // Convert semicolon-separated HPO terms to command line arguments
    def hpo_args = ''
    if (hpo_terms && hpo_terms != '') {
        hpo_args = hpo_terms.split(';').collect { "--phenotype-term ${it}" }.join(' ')
    }
    
    def output_formats = task.ext.output_format ?: 'html'
    
    """
    java -jar /app/svanna-cli.jar prioritize \\
        --data-directory ${data_directory} \\
        --output-format ${output_formats} \\
        --vcf ${vcf} \\
        ${hpo_args} \\
        --out-dir . \\
        --prefix ${prefix} \\
        --n-threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svanna: \$(java -jar /app/svanna-cli.jar --version 2>&1 | head -n1 | sed 's/.*v//' | sed 's/ .*//' || echo "1.0.4")
        java: \$(java -version 2>&1 | head -n1 | sed 's/.*version "//' | sed 's/".*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_formats = task.ext.output_format ?: 'html'
    
    """
    if [[ "${output_formats}" == *"html"* ]]; then
        touch ${prefix}.html
    fi
    if [[ "${output_formats}" == *"csv"* ]]; then
        touch ${prefix}.csv
    fi
    if [[ "${output_formats}" == *"vcf"* ]]; then
        touch ${prefix}.vcf.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svanna: "1.0.4"
        java: "17.0.12"
    END_VERSIONS
    """
}