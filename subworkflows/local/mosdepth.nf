//
// Subworkflow for running MOSDEPTH depth analysis
//

include { MOSDEPTH } from '../../modules/nf-core/mosdepth/main.nf'

// This workflow takes BAM files (with optional FASTA reference) and runs MOSDEPTH to calculate depth of coverage.

workflow mosdepth_subworkflow {
    take:
    bam_bai_bed // channel: [ val(meta), path(bam), path(bai), path(bed) ]
    fasta       // channel: [ val(meta2), path(fasta) ] - optional

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Run MOSDEPTH
    //
    MOSDEPTH (
        bam_bai_bed,
        fasta
    )
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())


    

    emit:
    global_txt     = MOSDEPTH.out.global_txt     // channel: [ val(meta), path(txt) ]
    summary_txt    = MOSDEPTH.out.summary_txt    // channel: [ val(meta), path(txt) ]
    regions_txt    = MOSDEPTH.out.regions_txt    // channel: [ val(meta), path(txt) ]
    per_base_d4    = MOSDEPTH.out.per_base_d4    // channel: [ val(meta), path(d4) ]
    per_base_bed   = MOSDEPTH.out.per_base_bed   // channel: [ val(meta), path(bed.gz) ]
    per_base_csi   = MOSDEPTH.out.per_base_csi   // channel: [ val(meta), path(csi) ]
    regions_bed    = MOSDEPTH.out.regions_bed    // channel: [ val(meta), path(bed.gz) ]
    regions_csi    = MOSDEPTH.out.regions_csi    // channel: [ val(meta), path(csi) ]
    thresholds_bed = MOSDEPTH.out.thresholds_bed // channel: [ val(meta), path(bed.gz) ]
    thresholds_csi = MOSDEPTH.out.thresholds_csi // channel: [ val(meta), path(csi) ]
    versions       = ch_versions                 // channel: [ versions.yml ]
}