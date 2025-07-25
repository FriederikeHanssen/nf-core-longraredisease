process TRUVARI_COLLAPSE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "biocontainers/truvari:5.3.0--pyhdfd78af_0"

    input:
    tuple val(meta), path(vcf), path(tbi)
    val(refdist) // Reference distance for merging
    val(pctsim)  // Percent similarity for merging
    val(pctseq)  // Percent sequence for merging
    val(no_consolidate) // Flag to disable consolidation
    val(passonly) // Flag to only keep PASS variants
    val(dup_to_ins) // Flag to treat Duplications as Insertions


    output:
    tuple val(meta), path("*.merged.vcf")         , emit: merged_vcf
    tuple val(meta), path("*.collapsed.vcf")      , emit: collapsed_vcf
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    

    """
    truvari collapse \\
        -p=0 \\
        --intra \\
        -i ${vcf} \\
        -o ${prefix}.merged.vcf \\
        -r ${refdist} \\
        -P ${pctsim} \\
        --pctseq ${pctseq} \\
        ${no_consolidate} \\
        ${passonly} \\
        ${dup_to_ins} \\
        -c ${prefix}.collapsed.vcf \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        truvari: \$(echo \$(truvari version 2>&1) | sed 's/^Truvari v//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged.vcf
    touch ${prefix}.collapsed.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        truvari: \$(echo \$(truvari version 2>&1) | sed 's/^Truvari v//' )
    END_VERSIONS
    """
}