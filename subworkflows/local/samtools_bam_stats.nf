include { BAM_STATS_SAMTOOLS } from '../../subworkflows/nf-core/bam_stats_samtools/main'

workflow bam_stats_subworkflow {

    take:
    ch_bam_bai
    ch_fasta

    main:

    // Run samtools stats
    BAM_STATS_SAMTOOLS(ch_bam_bai, ch_fasta)

    ch_stats     = BAM_STATS_SAMTOOLS.out.stats
    ch_flagstat  = BAM_STATS_SAMTOOLS.out.flagstat
    ch_idxstats  = BAM_STATS_SAMTOOLS.out.idxstats

    emit:
    stats      = ch_stats
    flagstat   = ch_flagstat
    idxstats   = ch_idxstats
}