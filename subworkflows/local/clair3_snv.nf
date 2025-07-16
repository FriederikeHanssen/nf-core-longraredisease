// This workflow is for clair3

include { CLAIR3 } from '../../modules/nf-core/clair3/main.nf'
include { FIX_CLAIR3} from '../../modules/local/fix_clair3_vcf/main'

workflow clair3_snv_subworkflow {
    take:
    input_bam    // channel: tuple(val(meta), path(bam), path(bai))
    fasta        // channel: tuple(val(meta2), path(fasta))
    fasta_fai    // channel: tuple(val(meta3), path(fai)) - optional

    main:
    ch_versions = Channel.empty()

    CLAIR3(
        input_bam,    // tuple(meta, bam, bai)
        fasta,        // tuple(meta2, fasta)
        fasta_fai     // tuple(meta3, fai)
    )

    // Join VCF and TBI to create the tuple format BCFTOOLS_VIEW expects
    

    FIX_CLAIR3(
        CLAIR3.out.vcf,    // path to VCF file
        CLAIR3.out.tbi,    // path to TBI file
    )       
    

    ch_versions = ch_versions.mix(CLAIR3.out.versions)

    emit:
    vcf        = FIX_CLAIR3.out.vcf
    tbi        = FIX_CLAIR3.out.tbi
    full_vcf   = CLAIR3.out.full_vcf
    full_tbi   = CLAIR3.out.full_tbi
    pileup_vcf = CLAIR3.out.pileup_vcf
    pileup_tbi = CLAIR3.out.pileup_tbi
    // phased_vcf = CLAIR3.out.phased_vcf
    // phased_tbi = CLAIR3.out.phased_tbi
    versions = ch_versions
}