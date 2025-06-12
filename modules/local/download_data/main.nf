process DOWNLOAD_DATA {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), val(dummy)  // meta map and dummy input to trigger process

    output:
    tuple val(meta), path("wf_str_repeat.bed"), emit: str_bed_file
    tuple val(meta), path("tandem_repeats.bed"), emit: tandem_repeats_bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Download the files using wget (available in most base systems)
    wget \\
        $args \\
        -O tandem_repeats.bed \\
        https://raw.githubusercontent.com/fritzsedlazeck/Sniffles/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed

    wget \\
        $args \\
        -O wf_str_repeat.bed \\
        https://raw.githubusercontent.com/epi2me-labs/wf-human-variation/master/data/wf_str_repeats.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -n1 | sed 's/GNU Wget //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch wf_str_repeat.bed
    touch tandem_repeats.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -n1 | sed 's/GNU Wget //')
    END_VERSIONS
    """
}