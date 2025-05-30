include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'

workflow bam_index_workflow {

    take:
    input_bam // channel: [ meta, bam, bai ]

    main:
    ch_versions = Channel.empty()

    // Run samtools index
    SAMTOOLS_INDEX (
        input_bam
    )

    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    emit:
    bai      = SAMTOOLS_INDEX.out.bai

    versions = ch_versions
}
