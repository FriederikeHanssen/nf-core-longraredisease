//
// Convert merged BAM files to extract unmapped/other reads using samtools fastq with methylation tags
//

include { SAMTOOLS_FASTQ } from '../../modules/nf-core/samtools/fastq/main'

workflow bam_to_fastq_workflow {
    
    take:
    ch_bam_files    // channel: [meta, bam] - output from bam_merge_workflow
    
    main:
    ch_versions = Channel.empty()
    
    // Convert BAM to extract other/unmapped reads with methylation tags
    // This will use: samtools fastq -T MM,ML,mv --threads 4 -0 output_other.fastq.gz input.bam
    SAMTOOLS_FASTQ(
        ch_bam_files,
        false  // interleave parameter set to false
    )
    
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions)
    
    emit:
    fastq       = SAMTOOLS_FASTQ.out.fastq       // channel: [meta, [fastq_1, fastq_2]] - paired-end files
    interleaved = SAMTOOLS_FASTQ.out.interleaved // channel: [meta, interleaved.fastq] - interleaved file
    singleton   = SAMTOOLS_FASTQ.out.singleton   // channel: [meta, singleton.fastq.gz] - singleton reads
    other       = SAMTOOLS_FASTQ.out.other       // channel: [meta, other.fastq.gz] - unmapped/other reads            // channel: [versions.yml]
}