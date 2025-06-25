include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_CUTESV } from '../../modules/nf-core/bedtools/intersect/main.nf'

workflow filterbyintersect_cutesv_subworkflow {

    take:
    ch_vcf_bed // channel: [ val(meta), path(intervals1) ]
    chrom_sizes // channel: [ val(meta3), path(chrom_sizes) ] - optional

    main:
    BEDTOOLS_INTERSECT_CUTESV (
        ch_vcf_bed,
        chrom_sizes
    )

    emit:
    intersect = BEDTOOLS_INTERSECT_CUTESV.out.intersect // channel: [ val(meta), path(intersect) ]
    versions  = BEDTOOLS_INTERSECT_CUTESV.out.versions  // channel: [ versions.yml ]
}