include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'

workflow samtools_index_subworkflow {

    take:
    ch_sorted_bam // channel: [meta, bam] - output from minimap2_align_workflow

    main:
    // Index BAM files using samtools index

    SAMTOOLS_INDEX(ch_sorted_bam)

    

    emit:
    bai   = SAMTOOLS_INDEX.out.bai // channel: [meta, bam.bai] - indexed BAM files
}