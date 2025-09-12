process HIFICNV {
    tag "$meta.id"
    label 'process_medium'

    container "quay.io/biocontainers/hificnv:1.0.0--h9ee0642_0"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(reference)
    path(exclude_bed)

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.depth.bw")  , emit: depth_bw
    tuple val(meta), path("*.bedgraph")  , emit: cnval
    tuple val(meta), path("*.log")       , emit: log
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def exclude = exclude_bed ? "--exclude ${exclude_bed}" : ""

    """
    hificnv \\
        $args \\
        --ref ${reference} \\
        --bam ${bam} \\
        $exclude \\
        --output-prefix ${prefix} \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hificnv: \$(hificnv -V | sed 's/hificnv //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.depth.bw
    touch ${prefix}.copynum.bedgraph
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hificnv: \$(hificnv -V | sed 's/hificnv //')
    END_VERSIONS
    """
}