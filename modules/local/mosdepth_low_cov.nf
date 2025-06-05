process MOSDEPTH_LOW_COV_FILTER {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(bed_gz)

    output:
    tuple val(meta), path("lowcov.bed"), emit: lowcov_bed

    script:
    """
    gunzip -c $bed_gz \
        | awk '\$4 < 10 { print \$1"\t"\$2"\t"\$3 }' > lowcov.bed
    """
}