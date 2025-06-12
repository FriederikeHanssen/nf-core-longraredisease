#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//
// Subworkflow to identify STR repeats in a genome using STAGLR
//

include { STRAGLR } from '../../modules/local/straglr/main'

workflow straglr_genotype_workflow {
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
    
    ch_versions = ch_versions.mix(STRAGLR.out.versions.first())
    
    emit:
    vcf      = STRAGLR.out.vcf      // channel: [ val(meta), path(vcf) ]
    versions = ch_versions         // channel: path(versions.yml)
}