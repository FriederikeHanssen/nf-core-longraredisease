//
// Subworkflow for running MOSDEPTH depth analysis
//

include { MOSDEPTH } from '../../modules/nf-core/mosdepth/main.nf'
include { MOSDEPTH_LOW_COV_FILTER} from '../../modules/local/mosdepth_low_cov/main.nf'
include { GUNZIP } from '../../modules/nf-core/gunzip/main.nf'
// This workflow takes BAM files (with optional FASTA reference) and runs MOSDEPTH to calculate depth of coverage.

workflow mosdepth_workflow {
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

    // Extract low coverage regions (<10x)
    MOSDEPTH_LOW_COV_FILTER(MOSDEPTH.out.per_base_bed)
    

    emit:
    global_txt     = MOSDEPTH.out.global_txt     // channel: [ val(meta), path(txt) ]
    summary_txt    = MOSDEPTH.out.summary_txt    // channel: [ val(meta), path(txt) ]
    regions_txt    = MOSDEPTH.out.regions_txt    // channel: [ val(meta), path(txt) ]
    per_base_d4    = MOSDEPTH.out.per_base_d4    // channel: [ val(meta), path(d4) ]
    per_base_bed   = MOSDEPTH.out.per_base_bed   // channel: [ val(meta), path(bed.gz) ]
    per_base_csi   = MOSDEPTH.out.per_base_csi   // channel: [ val(meta), path(csi) ]
    regions_bed    = MOSDEPTH.out.regions_bed    // channel: [ val(meta), path(bed.gz) ]
    regions_csi    = MOSDEPTH.out.regions_csi    // channel: [ val(meta), path(csi) ]
    quantized_bed  = MOSDEPTH.out.quantized_bed  // channel: [ val(meta), path(bed.gz) ]
    quantized_csi  = MOSDEPTH.out.quantized_csi  // channel: [ val(meta), path(csi) ]
    thresholds_bed = MOSDEPTH.out.thresholds_bed // channel: [ val(meta), path(bed.gz) ]
    thresholds_csi = MOSDEPTH.out.thresholds_csi // channel: [ val(meta), path(csi) ]
    lowcov_bed     = MOSDEPTH_LOW_COV_FILTER.out.lowcov_bed // channel: [ val(meta), path(bed) ]
    versions       = ch_versions                 // channel: [ versions.yml ]
}