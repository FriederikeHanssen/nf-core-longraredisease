// Subworkflow for running MOSDEPTH depth analysis

include { MOSDEPTH as MOSDEPTH_COV} from '../../modules/nf-core/mosdepth/main.nf'

// This workflow takes BAM files (with optional FASTA reference) and runs MOSDEPTH to calculate depth of coverage.

workflow mosdepth_cov_analysis_subworkflow {
    take:
    bam_bai_bed // channel: [ val(meta), path(bam), path(bai), path(bed) ]
    fasta       // channel: [ val(meta2), path(fasta) ] - optional

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Run MOSDEPTH
    //
    MOSDEPTH_COV (
        bam_bai_bed,
        fasta
    )
    ch_versions = ch_versions.mix(MOSDEPTH_COV.out.versions.first())


    

    emit:
    quantized_bed  = MOSDEPTH_COV.out.quantized_bed  // channel: [ val(meta), path(bed.gz) ]
    quantized_csi  = MOSDEPTH_COV.out.quantized_csi  // channel: [ val(meta), path(csi) ]
    versions       = ch_versions                 // channel: [ versions.yml ]
}