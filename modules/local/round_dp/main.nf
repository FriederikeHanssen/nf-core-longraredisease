process ROUND_DP {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.9"
    container "biocontainers/python:3.9--1"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Handle compressed input
    if [[ ${vcf} == *.gz ]]; then
        zcat ${vcf} | roundDP.py > ${prefix}.vcf
    else
        roundDP.py < ${vcf} > ${prefix}.vcf
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}