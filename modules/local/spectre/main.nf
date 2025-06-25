process SPECTRE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pip_ont-spectre:e37b5cd7d9231bc9' :
        'community.wave.seqera.io/library/pip_ont-spectre:b4773971f0688241' }"

    input:
    path(mosdepth_dir)
    val(bin_size)
    tuple val(meta), path(reference)
    path(snv_vcf)
    path(metadata_file)
    // path(blacklist)

    
    output:
    tuple val(meta), path("*.vcf")           , emit: vcf
    tuple val(meta), path("*.bed")           , emit: bed
    tuple val(meta), path("*.spc.gz")        , emit: spc
    tuple val(meta), path("*karyotype.txt")  , emit: karyo
    tuple val(meta), path("windows_stats")   , emit: winstats
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    spectre CNVCaller \\
        --coverage ${mosdepth_dir} \\
        --bin-size ${bin_size} \\
        --sample-id ${meta.id} \\
        --output-dir . \\
        --reference ${reference} \\
        --snv ${snv_vcf} \\
        --metadata ${metadata_file} \\
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