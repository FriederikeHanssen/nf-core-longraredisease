process SPECTRE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "wave.seqera.io/wt/7fb80fdd9784/wave/build:d9a1390fc3153e0e" // had to recreate the container as expired. Originally created using wave cli

    input:
    path(mosdepth_cov)
    tuple val(meta), path(reference)
    path(snv_vcf)
    path(metadata_file)
    path(blacklist)

    
    output:
    tuple val(meta), path("*.vcf.gz")        , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi")    , emit: index
    tuple val(meta), path("*.bed.gz")        , emit: bed
    tuple val(meta), path("*.bed.gz.tbi")    , emit: bed_index
    tuple val(meta), path("*.spc.gz")        , emit: spc
    // tuple val(meta), path("windows_stats")   , emit: winstats
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    spectre CNVCaller \\
        --coverage ${mosdepth_cov} \\
        --sample-id ${meta.id} \\
        --output-dir . \\
        --reference ${reference} \\
        --snv ${snv_vcf} \\
        --metadata ${metadata_file} \\
        --blacklist ${blacklist} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spectre: \$(spectre --version 2>&1 | grep -oP 'version \\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf
    touch ${prefix}.bed
    touch ${prefix}.spc.gz
    touch ${prefix}karyotype.txt
    touch versions.yml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spectre: \$(spectre --version 2>&1 | grep -oP 'version \\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """
}