process EXTRACT_BED {
    tag "$meta.id"
    label 'process_single'

    container 'community.wave.seqera.io/library/pip_awk:e0daab0638d06dfd'
    
    input:
    tuple val(meta), path(vcf_file)
    
    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"               , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk '!/^#/ {print \$1 "\\t" \$2 "\\t" \$3}' ${vcf_file} > ${prefix}.bed
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk --version | head -n1 | sed 's/GNU Awk //; s/,.*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk --version | head -n1 | sed 's/GNU Awk //; s/,.*//')
    END_VERSIONS
    """
}