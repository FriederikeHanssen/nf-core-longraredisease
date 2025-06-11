//
// Subworkflow for running SPECTRE CNV calling
//

include { SPECTRE } from '../../modules/local/spectre/main'

workflow spectre_cnv_workflow {
    take:
    ch_mosdepth_output    // channel: [ val(meta), path(mosdepth_dir) ]
    ch_reference          // channel: [ val(meta2), path(fasta) ]
    ch_snv_vcf            // channel: [ val(meta3), path(vcf) ]
    ch_metadata          // path to metadata file
    ch_blacklist         // path to blacklist file

    main:
    ch_versions = Channel.empty()
    //
    // MODULE: Run SPECTRE CNV calling
    //
    SPECTRE(
        ch_mosdepth_output,
        ch_reference,
        ch_snv_vcf,
        ch_metadata,
        ch_blacklist
    )
    ch_versions = ch_versions.mix(SPECTRE.out.versions.first())

    emit:
    vcf       = SPECTRE.out.vcf
    bed       = SPECTRE.out.bed
    spc       = SPECTRE.out.spc
    karyo     = SPECTRE.out.karyo
    winstats  = SPECTRE.out.winstats
    versions  = ch_versions 
}