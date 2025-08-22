// Convert merged BAM files to extract unmapped/other reads using samtools fastq with methylation tags

include { SAMTOOLS_MERGE } from '../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_FASTQ } from '../../modules/nf-core/samtools/fastq/main'

workflow bam2fastq_subworkflow {
    
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

    SAMTOOLS_FASTQ(
        SAMTOOLS_MERGE.out.bam,
        false  // interleave parameter set to false
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions)
    
    

emit:
    fastq       = SAMTOOLS_FASTQ.out.fastq       // channel: [meta, [fastq_1, fastq_2]] - paired-end files
    interleaved = SAMTOOLS_FASTQ.out.interleaved // channel: [meta, interleaved.fastq] - interleaved file
    singleton   = SAMTOOLS_FASTQ.out.singleton   // channel: [meta, singleton.fastq.gz] - singleton reads
    other       = SAMTOOLS_FASTQ.out.other       // channel: [meta, other.fastq.gz] - unmapped/other reads            // channel: [versions.yml]
    versions    = ch_versions                    // channel: [versions.yml] - versions of tools used in the workflow
}



