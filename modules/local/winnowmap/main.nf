process WINNOWMAP_ALIGN {
    tag "$meta.id"
    label 'process_high'

    container "community.wave.seqera.io/library/samtools_winnowmap:ed5f9fe3d48589a9"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(reference)
    path(kmers)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bam.bai"), emit: index
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def sort_args = task.ext.sort_args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def kmer_arg = kmers ? "-W ${kmers}" : ""
    def readgroup = "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ONT"
    
    """
    winnowmap \\
        ${kmer_arg} \\
        -R '${readgroup}' \\
        -t ${task.cpus} \\
        ${args} \\
        ${reference} \\
        ${reads} \\
        | samtools sort \\
            -@ ${task.cpus} \\
            ${sort_args} \\
            -o ${prefix}.bam \\
            -
    
    samtools index ${prefix}.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        winnowmap: \$(winnowmap --version 2>&1 | grep -o 'winnowmap-[0-9.]*' | sed 's/winnowmap-//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}