include { QDNASEQ } from '../../modules/local/qdnaseq/main.nf'

workflow cnv_qdnaseq_subworkflow {
    
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
    
    QDNASEQ(
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
    ch_versions = ch_versions.mix(QDNASEQ.out.versions)

    //
    // Collect all VCF files for downstream processing
    //
    ch_all_vcfs = QDNASEQ.out.calls_vcf
        .mix(QDNASEQ.out.segs_vcf)
        .mix(QDNASEQ.out.raw_calls_vcf)
        .mix(QDNASEQ.out.raw_segs_vcf)

    emit:
    calls_vcf     = QDNASEQ.out.calls_vcf     // channel: [meta, vcf]
    segs_vcf      = QDNASEQ.out.segs_vcf      // channel: [meta, vcf]
    raw_calls_vcf = QDNASEQ.out.raw_calls_vcf // channel: [meta, vcf]
    raw_segs_vcf  = QDNASEQ.out.raw_segs_vcf  // channel: [meta, vcf]
    plots         = QDNASEQ.out.plots         // channel: [meta, pdf]
    images        = QDNASEQ.out.images        // channel: [meta, png]
    tables        = QDNASEQ.out.tables        // channel: [meta, txt]
    bed           = QDNASEQ.out.bed           // channel: [meta, bed]
    seg           = QDNASEQ.out.seg           // channel: [meta, seg]
    all_vcfs      = ch_all_vcfs                  // channel: [meta, vcf]
    versions      = ch_versions                  // channel: versions.yml
}