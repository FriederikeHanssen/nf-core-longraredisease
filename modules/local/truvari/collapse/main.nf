process TRUVARI_COLLAPSE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "biocontainers/truvari:5.3.0--pyhdfd78af_0"

    input:
    tuple val(meta), path(vcf), path(tbi), path(bed)
    val(refdist)     // Reference distance for merging
    val(pctsim)      // Percent similarity for merging
    val(pctseq)      // Percent sequence for merging
    val(passonly)    // Boolean: Flag to only keep PASS variants
    val(dup_to_ins)  // Boolean: Flag to treat Duplications as Insertions

    output:
    tuple val(meta), path("*_merged.vcf")         , emit: merged_vcf
    tuple val(meta), path("*_collapsed.vcf")      , emit: collapsed_vcf
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // Convert boolean values to command-line flags
    def passonly_flag = passonly ? '--passonly' : ''
    def dup_to_ins_flag = dup_to_ins ? '--dup-to-ins' : ''
    def bed_flag = (bed && bed.name != 'NO_FILE') ? "--bed ${bed}" : ''

    """
    truvari collapse \\
        -i ${vcf} \\
        -o ${prefix}_merged.vcf \\
        -r ${refdist} \\
        -P ${pctsim} \\
        --pctseq ${pctseq} \\
        ${bed_flag} \\
        ${passonly_flag} \\
        ${dup_to_ins_flag} \\
        -c ${prefix}_collapsed.vcf \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        truvari: \$(echo \$(truvari version 2>&1) | sed 's/^Truvari v//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_merged.vcf
    touch ${prefix}_collapsed.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        truvari: \$(echo \$(truvari version 2>&1) | sed 's/^Truvari v//' )
    END_VERSIONS
    """
}