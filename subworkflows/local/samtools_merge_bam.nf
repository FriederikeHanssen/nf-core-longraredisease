//
// Subworkflow for merging BAM files
//

include { SAMTOOLS_MERGE } from '../../modules/nf-core/samtools/merge/main'

workflow samtools_merge_bam_subworkflow {
    
    take:
    ch_bam_files    // channel: [ meta, [ bam1, bam2, ... ] ]
    ch_fasta        // channel: [ meta, fasta ] (optional)
    ch_fai          // channel: [ meta, fai ] (optional)
    
    main:
    ch_versions = Channel.empty()
    
    // Handle optional inputs - use empty if not provided
    // ch_fasta_final = ch_fasta ?: ([[:], []])
    // ch_fai_final = ch_fai ?: ([[:], []])
    
    // Run samtools merge
    SAMTOOLS_MERGE (
        ch_bam_files,
        ch_fasta,
        ch_fai
    )
    
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)
    
    emit:
    bam      = SAMTOOLS_MERGE.out.bam
    versions = ch_versions
}
