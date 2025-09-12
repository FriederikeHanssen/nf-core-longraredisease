process HIFICNV {
    tag "$meta.id"
    label 'process_medium'

    container "community.wave.seqera.io/library/hificnv_htslib:ec2ecefb635052b8"

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

    # Fix sample name in VCF and rename files
    
    if [ -f "${prefix}.Sample0.vcf.gz" ]; then
    gunzip -c "${prefix}.Sample0.vcf.gz" | \\
    sed 's/Sample0/${meta.id}/g' | \\
    bgzip > "${prefix}.vcf.gz"
    
    # Index the new VCF
    
    tabix -p vcf "${prefix}.vcf.gz"
    
    # Remove original
    
    rm "${prefix}.Sample0.vcf.gz"
    fi
    
    # Rename other files
    
    if [ -f "${prefix}.Sample0.depth.bw" ]; then
    mv "${prefix}.Sample0.depth.bw" "${prefix}.depth.bw"
    fi
    
    if [ -f "${prefix}.Sample0.copynum.bedgraph" ]; then
    mv "${prefix}.Sample0.copynum.bedgraph" "${prefix}.copynum.bedgraph"
    fi
    
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