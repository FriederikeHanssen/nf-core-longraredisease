#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SVIM } from '../../modules/local/svim/main.nf'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_SVIM} from '../../modules/nf-core/tabix/bgziptabix/main.nf'

workflow svim_workflow {

    take:
    input         // tuple(val(meta), path(bam), path(bai))
    fasta         // tuple(val(meta), path(fasta))

    main:
    // 1. Call SVs with SVIM
    SVIM(input, fasta)

    // 2. Compress + index the VCF
    SVIM.out.vcf.map { meta, vcf -> tuple(meta, vcf) }
                .set { ch_vcf_to_index }

    TABIX_BGZIPTABIX_SVIM(ch_vcf_to_index)

    emit:
    vcf         = SVIM.out.vcf
    vcf_gz      = TABIX_BGZIPTABIX_SVIM.out.gz_tbi
    versions    = SVIM.out.versions
}
