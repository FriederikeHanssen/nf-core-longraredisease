// This workflow is for clair3

include { CLAIR3 } from '../../modules/nf-core/clair3/main.nf'
workflow clair3_snv_workflow {
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

    ch_versions = ch_versions.mix(CLAIR3.out.versions)

    emit:
    vcf      = CLAIR3.out.vcf
    tbi      = CLAIR3.out.tbi
    phased_vcf = CLAIR3.out.phased_vcf
    phased_tbi = CLAIR3.out.phased_tbi
    versions = ch_versions
}