// Subworkflow for running SPECTRE CNV calling

include { SPECTRE } from '../../modules/local/spectre/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_SPECTRE} from '../../modules/nf-core/tabix/bgziptabix/main.nf'

workflow spectre_cnv_subworkflow {
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
    index     = SPECTRE.out.index
    bed       = SPECTRE.out.bed
    bed_index = SPECTRE.out.bed_index
    spc       = SPECTRE.out.spc
    // winstats  = SPECTRE.out.winstats
    versions  = ch_versions 
}