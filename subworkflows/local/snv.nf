// This workflow is for clair3

include { CLAIR3 } from '../../modules/nf-core/clair3/main.nf'
include { CLAIR3_FIX} from '../../modules/local/clair3_fixVCF/main'
include { DEEPVARIANT_RUNDEEPVARIANT } from '../../modules/nf-core/deepvariant/rundeepvariant/main.nf'
include { BCFTOOLS_VIEW as BCFTOOLS_FILTER_CLAIR3 } from '../../modules/nf-core/bcftools/view/main.nf'
include { BCFTOOLS_VIEW as BCFTOOLS_FILTER_DEEPVARIANT } from '../../modules/nf-core/bcftools/view/main.nf'

workflow snv_subworkflow {
    take:
    ch_input_clair3         // channel: tuple(val(meta), path(bam), path(bai))
    fasta                   // channel: tuple(val(meta2), path(fasta))
    fai                     // channel: tuple(val(meta3), path(fai)) - optional
    deepvariant             // boolean, if true, will run deepvariant on the input bam file
    ch_input_deepvariant    // channel: [meta, bam, bai]
    filter_pass             // boolean, if true, will filter for PASS variants
    
    main:
    ch_versions = Channel.empty()

    // Run CLAIR3
    CLAIR3(
        ch_input_clair3,    // tuple(meta, bam, bai)
        fasta,              // tuple(meta2, fasta)
        fai                 // tuple(meta3, fai)
    )

    // Fix CLAIR3 VCF
    CLAIR3_FIX(
        CLAIR3.out.vcf,     // path to VCF file
        CLAIR3.out.tbi      // path to TBI file
    )

    // Collect versions
    ch_versions = ch_versions.mix(CLAIR3.out.versions)
    ch_versions = ch_versions.mix(CLAIR3_FIX.out.versions)

    // Handle CLAIR3 filtering
    if (filter_pass) {
        ch_clair3_vcf = CLAIR3_FIX.out.vcf
            .join(CLAIR3_FIX.out.tbi, by: 0)
         
        BCFTOOLS_FILTER_CLAIR3(
            ch_clair3_vcf,      // tuple(meta, vcf, tbi)
            Channel.value([]),   // empty channel for samples
            Channel.value([]),   // empty channel for regions
            Channel.value([])    // empty channel for filters
        )

        ch_final_clair3_vcf = BCFTOOLS_FILTER_CLAIR3.out.vcf
        ch_final_clair3_tbi = BCFTOOLS_FILTER_CLAIR3.out.tbi
        ch_versions = ch_versions.mix(BCFTOOLS_FILTER_CLAIR3.out.versions)
        
    } else {
        ch_final_clair3_vcf = CLAIR3_FIX.out.vcf
        ch_final_clair3_tbi = CLAIR3_FIX.out.tbi
    }

    // Handle DeepVariant
    if (deepvariant) {
        DEEPVARIANT_RUNDEEPVARIANT(
            ch_input_deepvariant,
            fasta,
            fai,
            [[:], []],
            [[:], []]
        )

        ch_deepvariant_versions = DEEPVARIANT_RUNDEEPVARIANT.out.versions

        if (filter_pass) {
            ch_deepvariant_vcf = DEEPVARIANT_RUNDEEPVARIANT.out.vcf
                .join(DEEPVARIANT_RUNDEEPVARIANT.out.vcf_tbi, by: 0)
            
            BCFTOOLS_FILTER_DEEPVARIANT(
                ch_deepvariant_vcf,
                Channel.value([]),
                Channel.value([]),  
                Channel.value([])
            )

            ch_final_deepvariant_vcf = BCFTOOLS_FILTER_DEEPVARIANT.out.vcf
            ch_final_deepvariant_tbi = BCFTOOLS_FILTER_DEEPVARIANT.out.tbi
            ch_versions = ch_versions.mix(BCFTOOLS_FILTER_DEEPVARIANT.out.versions)
            
        } else {
            ch_final_deepvariant_vcf = DEEPVARIANT_RUNDEEPVARIANT.out.vcf
            ch_final_deepvariant_tbi = DEEPVARIANT_RUNDEEPVARIANT.out.vcf_tbi
        }
        
    } else {
        // Create empty channels when DeepVariant is not run
        ch_final_deepvariant_vcf = Channel.empty()
        ch_final_deepvariant_tbi = Channel.empty()
        ch_deepvariant_versions = Channel.empty()
    }

    emit:
    clair3_vcf           = ch_final_clair3_vcf      // Filtered or unfiltered based on filter_pass
    clair3_tbi           = ch_final_clair3_tbi      // Corresponding index
    clair3_full_vcf      = CLAIR3.out.full_vcf      // Original full VCF
    clair3_full_tbi      = CLAIR3.out.full_tbi      // Original full TBI
    clair3_pileup_vcf    = CLAIR3.out.pileup_vcf    // Pileup VCF
    clair3_pileup_tbi    = CLAIR3.out.pileup_tbi    // Pileup TBI
    deepvariant_vcf      = ch_final_deepvariant_vcf // Filtered or unfiltered DeepVariant VCF
    deepvariant_tbi      = ch_final_deepvariant_tbi // Corresponding index
    clair3_versions      = ch_versions              // All collected versions
    deepvariant_versions = ch_deepvariant_versions  // DeepVariant versions
}