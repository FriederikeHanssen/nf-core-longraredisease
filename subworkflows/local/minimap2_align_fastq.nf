include { MINIMAP2_INDEX as MINIMAP2_INDEX_FASTQ } from '../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_FASTQ } from '../../modules/nf-core/minimap2/align/main'

workflow minimap2_align_fastq_subworkflow {

    take:
    ch_fasta        // channel: [meta, fasta]
    ch_fastq        // channel: [meta, [fastq files]]

    main:
    ch_versions = Channel.empty()

    // 1. Generate minimap2 index
    MINIMAP2_INDEX_FASTQ(ch_fasta)
    ch_versions = ch_versions.mix(MINIMAP2_INDEX_FASTQ.out.versions)

    // 2. Prepare channels for alignment
    // Create all combinations of fastq samples with the reference
    ch_fastq
        .combine(MINIMAP2_INDEX_FASTQ.out.index)
        .map { meta_sample, fastq, meta_ref, index ->
            [
                meta_sample,           // sample metadata
                fastq,                 // fastq files
                meta_ref,              // reference metadata  
                index                  // minimap2 index
            ]
        }
        .set { ch_align_input }

    // 3. Run minimap2 alignment
    // Set parameters for ONT long-read alignment
    bam_format = true
    bam_index_extension = 'bai'
    cigar_paf_format = false
    cigar_bam = false

    MINIMAP2_ALIGN_FASTQ(
        ch_align_input.map { meta_sample, fastq, meta_ref, index -> 
            [meta_sample, fastq] 
        },
        ch_align_input.map { meta_sample, fastq, meta_ref, index -> 
            [meta_ref, index] 
        },
        bam_format,
        bam_index_extension,
        cigar_paf_format,
        cigar_bam
    )
    
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN_FASTQ.out.versions)

    emit:
    bam      = MINIMAP2_ALIGN_FASTQ.out.bam     // channel: [meta, bam]
    index    = MINIMAP2_ALIGN_FASTQ.out.index   // channel: [meta, bai]
    versions = ch_versions                 // channel: [versions.yml]
}