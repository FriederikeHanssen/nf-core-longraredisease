include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_SNIFFLES } from '../../modules/nf-core/bedtools/intersect/main.nf'

workflow filterbyintersect_sniffles_subworkflow {

    take:
    ch_vcf_bed // channel: [ val(meta), path(intervals1) ]
    chrom_sizes // channel: [ val(meta3), path(chrom_sizes) ] - optional

    main:
    BEDTOOLS_INTERSECT_SNIFFLES (
        ch_vcf_bed,
        chrom_sizes
    )

    emit:
    intersect = BEDTOOLS_INTERSECT_SNIFFLES.out.intersect // channel: [ val(meta), path(intersect) ]
    versions  = BEDTOOLS_INTERSECT_SNIFFLES.out.versions  // channel: [ versions.yml ]
}