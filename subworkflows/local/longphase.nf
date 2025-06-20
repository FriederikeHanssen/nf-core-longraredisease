

include { LONGPHASE_PHASE } from '../../modules/nf-core/longphase/phase/main.nf'

workflow longphase_subworkflow {
    take:
    ch_bam // channel: tuple(val(meta), path(bam), path(bai))
    ch_fasta // channel: tuple(val(meta2), path(fasta))
    ch_fai   // channel: tuple(val(meta3), path(fai))

    main:
    ch_versions = Channel.empty()
    LONGPHASE_PHASE(
        ch_bam,  // tuple(meta, bam, bai)
        ch_fasta, // tuple(meta2, fasta)
        ch_fai    // tuple(meta3, fai)
    )

    ch_versions = ch_versions.mix(LONGPHASE_PHASE.out.versions)
    emit:
    vcf = LONGPHASE_PHASE.out.vcf
    sv_vcf= LONGPHASE_PHASE.out.sv_vcf
    versions = ch_versions
}