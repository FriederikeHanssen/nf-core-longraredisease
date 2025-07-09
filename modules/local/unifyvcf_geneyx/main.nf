process UNIFYVCF {
    tag "$meta.id"
    label 'process_low'
    
    conda "${moduleDir}/environment.yml"
    container "wave.seqera.io/wt/211639257add/wave/build:python-3.9_htslib--9c2949194826ffc2"

    input:
    tuple val(meta), path(sv_vcf)       // Single SV VCF file
    tuple val(meta1), path(cnv_vcf)     // Single CNV VCF file
    tuple val(meta2), path(repeat_vcf)  // Single STR VCF file
    val(modify_repeats)                 // Boolean: whether to modify repeat calls

    output:
    tuple val(meta), path("${prefix}_unified.vcf.gz"), emit: unified_vcf
    tuple val(meta), path("${prefix}_unified.vcf.gz.tbi"), emit: unified_tbi
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    
    // Handle optional VCF inputs - check if files exist and are not empty
    def sv_arg = sv_vcf  ? "-s ${sv_vcf}" : ''
    def cnv_arg = cnv_vcf ? "-c ${cnv_vcf}" : ''
    def repeat_arg = repeat_vcf ? "-r ${repeat_vcf}" : ''
    def modify_repeats_arg = modify_repeats ? '-modify' : ''
    
    """
    ONTUnifyVcf.py \\
        -o ${prefix}_unified.vcf \\
        ${sv_arg} \\
        ${cnv_arg} \\
        ${repeat_arg} \\
        ${modify_repeats_arg} \\
        ${args}

    tabix ${prefix}_unified.vcf.gz

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