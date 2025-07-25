// Subworkflow for running MOSDEPTH depth analysis

include { MOSDEPTH } from '../../modules/nf-core/mosdepth/main.nf'
include { MOSDEPTH as MOSDEPTH_COV} from '../../modules/nf-core/mosdepth/main.nf'

// This workflow takes BAM files (with optional FASTA reference) and runs MOSDEPTH to calculate depth of coverage.

workflow mosdepth_subworkflow {
    take:
    bam_bai_bed // channel: [ val(meta), path(bam), path(bai), path(bed) ]
    fasta       // channel: [ val(meta2), path(fasta) ] - optional

    
    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Run MOSDEPTH (always)
    //
    MOSDEPTH (
        bam_bai_bed,
        fasta
    )

    //
    // MODULE: Run MOSDEPTH_COV (always, but only use output if needed)
    //
    MOSDEPTH_COV (
        bam_bai_bed,
        fasta
    )

    ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())
    ch_versions = ch_versions.mix(MOSDEPTH_COV.out.versions.first())

    emit:
    global_txt     = MOSDEPTH.out.global_txt     
    summary_txt    = MOSDEPTH.out.summary_txt    
    regions_txt    = MOSDEPTH.out.regions_txt    
    per_base_d4    = MOSDEPTH.out.per_base_d4    
    per_base_bed   = MOSDEPTH.out.per_base_bed   
    per_base_csi   = MOSDEPTH.out.per_base_csi   
    regions_bed    = MOSDEPTH.out.regions_bed    
    regions_csi    = MOSDEPTH.out.regions_csi    
    thresholds_bed = MOSDEPTH.out.thresholds_bed 
    thresholds_csi = MOSDEPTH.out.thresholds_csi 
    quantized_bed  = MOSDEPTH_COV.out.quantized_bed  
    quantized_csi  = MOSDEPTH_COV.out.quantized_csi  
    versions       = ch_versions                 
}
