// Subworkflow to identify STR repeats in a genome using STAGLR

include { STRAGLR } from '../../modules/local/straglr/main'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_STRAGLR } from '../../modules/nf-core/bcftools/sort/main.nf'

workflow str_subworkflow {
    take:
    ch_bam_bai    // channel: [ val(meta), path(bam), path(bai) ]
    ch_reference  // channel: [ val(meta2), path(reference) ]
    ch_bed_file   // channel: path(bed_file)

    main:
    ch_versions = Channel.empty()
    
    STRAGLR(
        ch_bam_bai,
        ch_reference,
        ch_bed_file
    )
    
    BCFTOOLS_SORT_STRAGLR(
        STRAGLR.out.vcf
    )
    ch_versions = ch_versions.mix(STRAGLR.out.versions.first())
    
    emit:
    vcf      = BCFTOOLS_SORT_STRAGLR.out.vcf      // channel: [ val(meta), path(vcf) ]
    versions = ch_versions         // channel: path(versions.yml)
}