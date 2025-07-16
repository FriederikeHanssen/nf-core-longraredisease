// Copy Number Variant Detection Subworkflow


include { RUNQDNASEQ } from '../../modules/local/qdnaseq/main.nf'

workflow qdnaseq_cnv_subworkflow {
    
    take:
    ch_bam_bai      // channel: [meta, bam, bai]
    val_genome      // val: genome build (hg38/hg19)
    val_bin_size    // val: bin size in kbp
    val_method      // val: analysis method (cutoff/CGHcall)
    val_cutoff      // val: CN cutoff modifier
    val_cutoff_del  // val: deletion threshold
    val_cutoff_loss // val: loss threshold
    val_cutoff_gain // val: gain threshold
    val_cellularity // val: cellularity value

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Run QDNAseq copy number analysis
    //
    RUNQDNASEQ(
        ch_bam_bai,
        val_genome,
        val_bin_size,
        val_method,
        val_cutoff,
        val_cutoff_del,
        val_cutoff_loss,
        val_cutoff_gain,
        val_cellularity
    )
    ch_versions = ch_versions.mix(RUNQDNASEQ.out.versions)

    //
    // Collect all VCF files for downstream processing
    //
    ch_all_vcfs = RUNQDNASEQ.out.calls_vcf
        .mix(RUNQDNASEQ.out.segs_vcf)
        .mix(RUNQDNASEQ.out.raw_calls_vcf)
        .mix(RUNQDNASEQ.out.raw_segs_vcf)

    emit:
    calls_vcf     = RUNQDNASEQ.out.calls_vcf     // channel: [meta, vcf]
    segs_vcf      = RUNQDNASEQ.out.segs_vcf      // channel: [meta, vcf]
    raw_calls_vcf = RUNQDNASEQ.out.raw_calls_vcf // channel: [meta, vcf]
    raw_segs_vcf  = RUNQDNASEQ.out.raw_segs_vcf  // channel: [meta, vcf]
    plots         = RUNQDNASEQ.out.plots         // channel: [meta, pdf]
    images        = RUNQDNASEQ.out.images        // channel: [meta, png]
    tables        = RUNQDNASEQ.out.tables        // channel: [meta, txt]
    bed           = RUNQDNASEQ.out.bed           // channel: [meta, bed]
    seg           = RUNQDNASEQ.out.seg           // channel: [meta, seg]
    all_vcfs      = ch_all_vcfs                  // channel: [meta, vcf]
    versions      = ch_versions                  // channel: versions.yml
}