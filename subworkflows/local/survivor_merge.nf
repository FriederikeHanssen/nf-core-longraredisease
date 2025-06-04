#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SURVIVOR_MERGE } from '../../modules/nf-core/survivor/merge/main.nf'

workflow survivor_merge_workflow {

    take:
    vcfs         // Channel: tuple(meta, List[VCF file])
    max_dist     // int
    min_callers  // int
    type_flag    // int
    strand_flag  // int
    est_dist     // int
    min_size     // int

    main:
    vcfs
        .map { meta, files ->
            tuple(meta, files)
        }
        .set { ch_vcf_input }

    SURVIVOR_MERGE(
        ch_vcf_input,
        max_dist,
        min_callers,
        type_flag,
        strand_flag,
        est_dist,
        min_size
    )

    emit:
    vcf      = SURVIVOR_MERGE.out.vcf
    versions = SURVIVOR_MERGE.out.versions
}
