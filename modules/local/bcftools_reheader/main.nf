process BCFTOOLS_REHEADER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.{vcf,vcf.gz,bcf,bcf.gz}"), emit: vcf
    tuple val(meta), path("*.{csi,tbi}")              , emit: index, optional: true
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = vcf.name.endsWith('.gz') ? 'vcf.gz' : 'vcf'

    """
    bcftools \\
        reheader \\
        -s <(echo -e "${meta.id}") \\
        $args \\
        $vcf \\
        -o ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = vcf.name.endsWith('.gz') ? 'vcf.gz' : 'vcf'
    
    """
    touch ${prefix}.${extension}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}