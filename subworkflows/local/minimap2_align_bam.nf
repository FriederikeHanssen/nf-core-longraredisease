//Alignment with MINIMAP2

include { MINIMAP2_INDEX as MINIMAP2_INDEX_BAM } from '../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_BAM } from '../../modules/nf-core/minimap2/align/main'

workflow minimap2_align_bam_subworkflow {

    take:
    ch_fasta
    ch_fastq

    main:

    // 1. Generate index (runs once)
    MINIMAP2_INDEX_BAM(ch_fasta)
    ch_minimap_index = MINIMAP2_INDEX_BAM.out.index
    minimap2_version = MINIMAP2_INDEX_BAM.out.versions

    // 2. Update meta so it is updated with sample id
    ch_align_input = ch_fastq
        .combine(ch_minimap_index)
        .map { meta_sample, reads, meta_ref, index -> 
            // Create tuple with sample meta for both reads and index
            [meta_sample, reads, meta_sample, index]
        }

    // 3. Run alignment
    def bam_format = true
    def bam_index_extension = 'bai'
    def cigar_paf_format = false
    def cigar_bam = false

    MINIMAP2_ALIGN_BAM(
        ch_align_input.map { meta_sample, reads, meta_index, index -> 
            [meta_sample, reads] 
        },
        ch_align_input.map { meta_sample, reads, meta_index, index -> 
            [meta_index, index]  // Use sample meta for index too
        },
        bam_format,
        bam_index_extension,
        cigar_paf_format,
        cigar_bam
    )

    ch_sorted_bam = MINIMAP2_ALIGN_BAM.out.bam
    ch_sorted_bai = MINIMAP2_ALIGN_BAM.out.index

    emit:
    ch_sorted_bam
    ch_sorted_bai
    versions = minimap2_version
}