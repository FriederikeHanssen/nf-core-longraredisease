// Alignment with MINIMAP2 or WINNOWMAP
include { MINIMAP2_INDEX  } from '../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'
include { WINNOWMAP_ALIGN } from '../../modules/local/winnowmap/main'

workflow alignment_subworkflow {

    take:
    ch_fasta
    ch_fastq
    ch_kmers       // Optional: winnowmap k-mer file (empty channel if not needed)

    main:

    ch_versions = Channel.empty()
    
    if (params.use_winnowmap) {
        // WINNOWMAP alignment
        WINNOWMAP_ALIGN(
            ch_fastq,
            ch_fasta,
            ch_kmers
        )
        
        ch_sorted_bam = WINNOWMAP_ALIGN.out.bam
        ch_sorted_bai = WINNOWMAP_ALIGN.out.index
        ch_versions = ch_versions.mix(WINNOWMAP_ALIGN.out.versions)
        
    } else {
        // MINIMAP2 alignment (default)
        
        // 1. Generate index (runs once)
        MINIMAP2_INDEX(ch_fasta)
        ch_minimap_index = MINIMAP2_INDEX.out.index
        ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)

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

        MINIMAP2_ALIGN(
            ch_align_input.map { meta_sample, reads, meta_index, index -> 
                [meta_sample, reads] 
            },
            ch_align_input.map { meta_sample, reads, meta_index, index -> 
                [meta_index, index]  // Use sample meta for index too as otherwise runs one sample at a time 
            },
            bam_format,
            bam_index_extension,
            cigar_paf_format,
            cigar_bam
        )

        ch_sorted_bam = MINIMAP2_ALIGN.out.bam
        ch_sorted_bai = MINIMAP2_ALIGN.out.index
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
    }

    emit:
    bam = ch_sorted_bam
    bai = ch_sorted_bai
    versions = ch_versions
}