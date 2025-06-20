#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SURVIVOR_MERGE } from '../../modules/nf-core/survivor/merge/main.nf'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_SURVIVOR} from '../../modules/nf-core/tabix/bgziptabix/main.nf'

workflow survivor_merge_sv_subworkflow {

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

// max_dist = 1000
    SURVIVOR_MERGE(
        ch_vcf_input,
        max_dist,
        min_callers,
        type_flag,
        strand_flag,
        est_dist,
        min_size
    )

    SURVIVOR_MERGE.out.vcf.map { meta, vcf -> tuple(meta, vcf) }
                  .set { ch_vcf_to_index }

    // TABIX_BGZIPTABIX_SURVIVOR(ch_vcf_to_index)

    emit:
    vcf      = SURVIVOR_MERGE.out.vcf
    // vcf_gz   = TABIX_BGZIPTABIX_SURVIVOR.out.gz_tbi
    versions = SURVIVOR_MERGE.out.versions
}
