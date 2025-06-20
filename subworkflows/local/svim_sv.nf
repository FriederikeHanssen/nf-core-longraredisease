#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SVIM } from '../../modules/local/svim/main.nf'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_SVIM } from '../../modules/nf-core/bcftools/sort/main.nf'
// include { TABIX_TABIX as TABIX_TABIX_SVIM} from '../../modules/nf-core/tabix/tabix/main.nf'

workflow svim_sv_subworkflow {

    take:
    input         // tuple(val(meta), path(bam), path(bai))
    fasta         // tuple(val(meta), path(fasta))

    main:
    // 1. Call SVs with SVIM
    SVIM(input, fasta)

    // 2. Compress + index the VCF
    SVIM.out.vcf.map { meta, vcf -> tuple(meta, vcf) }
                .set { ch_vcf_to_index }

   BCFTOOLS_SORT_SVIM(ch_vcf_to_index)

//    TABIX_TABIX_SVIM(BCFTOOLS_SORT_SVIM.out.vcf)

//     TABIX_BGZIPTABIX_SVIM(BCFTOOLS_SORT.out.vcf)

    emit:
    vcf         = SVIM.out.vcf
    vcf_gz      = BCFTOOLS_SORT_SVIM.out.vcf
    // tbi          = TABIX_TABIX_SVIM.out.tbi
    versions    = SVIM.out.versions
}
