process SVIM {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svim:2.0.0--pyhdfd78af_0' :
        'biocontainers/svim:2.0.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def outdir = "${meta.id}_svim"

    """
    svim alignment \
    ${outdir} \\
    ${bam} \\
    ${fasta} \\
    $args

    mv ${outdir}/variants.vcf ${meta.id}.vcf

    cat <<-END_VERSIONS > versions.yml 
    "${task.process}":
        svim: \$(svim --version 2>&1 | sed 's/svim //g')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.vcf"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svim: \$(svim --version 2>&1 | sed 's/svim //g')
    END_VERSIONS
    """
}
