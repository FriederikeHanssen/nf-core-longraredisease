//
// Subworkflow for running SPECTRE CNV calling
//

include { SPECTRE } from '../../modules/local/spectre/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_SPECTRE} from '../../modules/nf-core/tabix/bgziptabix/main.nf'

workflow spectre_cnv_subworkflow {
    take:
    ch_mosdepth_output    // channel: [ val(meta), path(mosdepth_dir) ]
    bin_size              // must be the same as mosdepth
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
        bin_size,
        ch_reference,
        ch_snv_vcf,
        ch_metadata,
        ch_blacklist
    )
    ch_versions = ch_versions.mix(SPECTRE.out.versions.first())
    SPECTRE.out.vcf.map { meta, vcf -> tuple(meta, vcf) }
                .set { ch_vcf_to_index }

    TABIX_BGZIPTABIX_SPECTRE(ch_vcf_to_index)


    emit:
    vcf       = SPECTRE.out.vcf
    vcf_gz     = TABIX_BGZIPTABIX_SPECTRE.out.gz_tbi
    bed       = SPECTRE.out.bed
    spc       = SPECTRE.out.spc
    karyo     = SPECTRE.out.karyo
    winstats  = SPECTRE.out.winstats
    versions  = ch_versions 
}