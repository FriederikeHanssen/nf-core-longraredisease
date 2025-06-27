process UNIFY_VCFS {
    tag "$meta.id"
    label 'process_low'
    
    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(sv_vcfs)      // Multiple SV VCF files as a list
    tuple val(meta2), path(cnv_vcf)     // Single CNV VCF file (optional)
    tuple val(meta3), path(repeat_vcf)  // Single repeat VCF file (optional)

    output:
    tuple val(meta), path("${prefix}_unified.vcf"), emit: unified_vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    
    // Build the command arguments for multiple SV files
    def sv_args = ""
    if (sv_vcfs instanceof List && sv_vcfs.size() > 0) {
        // Multiple SV files
        sv_args = sv_vcfs.findAll { it.name != 'OPTIONAL_FILE' }
                         .collect { "--svPaths $it" }
                         .join(' ')
    } else if (sv_vcfs && sv_vcfs.name != 'OPTIONAL_FILE') {
        // Single SV file
        sv_args = "--svPaths $sv_vcfs"
    }
    
    def cnv_arg = (cnv_vcf && cnv_vcf.name != 'OPTIONAL_FILE') ? "--cnvPath $cnv_vcf" : ''
    def repeat_arg = (repeat_vcf && repeat_vcf.name != 'OPTIONAL_FILE') ? "--repeatPath $repeat_vcf" : ''
    
    """
    unify_vcf.py \\
        --outputPath ${prefix}_unified.vcf \\
        $sv_args \\
        $cnv_arg \\
        $repeat_arg \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_unified.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}