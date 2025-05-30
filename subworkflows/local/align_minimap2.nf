/*
 * Alignment with MINIMAP2
 */

include { MINIMAP2_INDEX } from '../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../../subworkflows/nf-core/bam_stats_samtools/main'

workflow ALIGN_MINIMAP2 {

    take:
    ch_fasta
    ch_fastq

    main:

    // 1. Generate index
    MINIMAP2_INDEX ( ch_fasta )
    ch_minimap_index   = MINIMAP2_INDEX.out.index
    minimap2_version   = MINIMAP2_INDEX.out.versions

    // 2. Map reads to reference
    ch_fastq
        .map { meta, reads -> tuple(meta, reads) }
        .set { ch_reads }

    ch_minimap_index
        .map { meta, ref -> tuple(meta, ref) }
        .set { ch_ref }

    // These can also be set dynamically if needed
    def bam_format = true
    def bam_index_extension = 'bai'
    def cigar_paf_format = false
    def cigar_bam = false

    MINIMAP2_ALIGN(
        ch_reads,
        ch_ref,
        bam_format,
        bam_index_extension,
        cigar_paf_format,
        cigar_bam
    )

    ch_sorted_bam     = MINIMAP2_ALIGN.out.bam
    ch_sorted_bai     = MINIMAP2_ALIGN.out.index
    minimap2_version  = MINIMAP2_ALIGN.out.versions

    // 3. Index BAM
    SAMTOOLS_INDEX(ch_sorted_bam)
    ch_sorted_bai = SAMTOOLS_INDEX.out.bai
    samtools_version = SAMTOOLS_INDEX.out.versions

    ch_sorted_bam
        .join(ch_sorted_bai, by: 0)
        .set { ch_bam_bai }

    // 4. Run samtools stats
    BAM_STATS_SAMTOOLS(ch_bam_bai, ch_fasta)
    ch_stats     = BAM_STATS_SAMTOOLS.out.stats
    ch_flagstat  = BAM_STATS_SAMTOOLS.out.flagstat
    ch_idxstats  = BAM_STATS_SAMTOOLS.out.idxstats

    emit:
    ch_minimap_index
    minimap2_version
    ch_sorted_bam
    ch_sorted_bai
    samtools_version
    ch_stats
    ch_flagstat
    ch_idxstats
}
