//
// Run DeepVariant variant calling
//

include { DEEPVARIANT_RUNDEEPVARIANT } from '../../modules/nf-core/deepvariant/rundeepvariant/main.nf'

workflow deepvariant_snv_subworkflow {
    take:
    ch_input    // channel: [meta, bam, bai]
    ch_fasta    // channel: [meta, fasta]
    ch_fai      // channel: [meta, fai]
    ch_gzi      // channel: [meta, gzi]
    ch_par_bed  // channel: [meta, par_bed]

    main:
    ch_versions = Channel.empty()



    // Run DeepVariant
    DEEPVARIANT_RUNDEEPVARIANT(
        ch_input,
        ch_fasta,
        ch_fai,
        ch_gzi,
        ch_par_bed
    )

    ch_versions = ch_versions.mix(DEEPVARIANT_RUNDEEPVARIANT.out.versions)

    emit:
    vcf      = DEEPVARIANT_RUNDEEPVARIANT.out.vcf      // channel: [meta, vcf]
    vcf_tbi  = DEEPVARIANT_RUNDEEPVARIANT.out.vcf_tbi  // channel: [meta, vcf_tbi]
    gvcf     = DEEPVARIANT_RUNDEEPVARIANT.out.gvcf     // channel: [meta, gvcf]
    gvcf_tbi = DEEPVARIANT_RUNDEEPVARIANT.out.gvcf_tbi // channel: [meta, gvcf_tbi]
    versions = ch_versions                              // channel: [versions.yml]
}