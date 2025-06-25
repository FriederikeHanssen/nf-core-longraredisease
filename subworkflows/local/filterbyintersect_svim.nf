include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_SVIM } from '../../modules/nf-core/bedtools/intersect/main.nf'

workflow filterbyintersect_svim_subworkflow {

    take:
    ch_vcf_bed // channel: [ val(meta), path(intervals1) ]
    chrom_sizes // channel: [ val(meta3), path(chrom_sizes) ] - optional

    main:
    BEDTOOLS_INTERSECT_SVIM (
        ch_vcf_bed,
        chrom_sizes
    )

    emit:
    intersect = BEDTOOLS_INTERSECT_SVIM.out.intersect // channel: [ val(meta), path(intersect) ]
    versions  = BEDTOOLS_INTERSECT_SVIM.out.versions  // channel: [ versions.yml ]
}