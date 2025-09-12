// Copy Number Variant Detection Subworkflow

include { SPECTRE } from '../../modules/local/spectre/main'
include { ROUND_DP } from '../../modules/local/round_dp/main'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_SPECTRE } from '../../modules/nf-core/bcftools/sort'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_SPECTRE} from '../../modules/nf-core/tabix/bgziptabix/main.nf'


workflow cnv_spectre_subworkflow {
    take:
    ch_mosdepth_output    // channel: [ val(meta), path(mosdepth_dir) ]
    ch_reference          // channel: [ val(meta2), path(fasta) ]
    ch_snv_vcf            // channel: [ val(meta3), path(vcf) ]
    ch_metadata          // path to metadata file
    ch_blacklist         // path to blacklist file

    main:

    //
    // MODULE: Run SPECTRE CNV calling
    //
    ch_versions = Channel.empty()

    SPECTRE(
        ch_mosdepth_output,
        ch_reference,
        ch_snv_vcf,
        ch_metadata,
        ch_blacklist
    )

    ch_versions = ch_versions.mix(SPECTRE.out.versions)
    
    ROUND_DP(SPECTRE.out.vcf)

    ch_versions = ch_versions.mix(ROUND_DP.out.versions)

    BCFTOOLS_SORT_SPECTRE(ROUND_DP.out.vcf)
    
   ch_versions = ch_versions.mix(BCFTOOLS_SORT_SPECTRE.out.versions)

    emit:
    vcf       = BCFTOOLS_SORT_SPECTRE.out.vcf
    tbi       = BCFTOOLS_SORT_SPECTRE.out.tbi
    bed       = SPECTRE.out.bed
    bed_index = SPECTRE.out.bed_index
    spc       = SPECTRE.out.spc
    winstats  = SPECTRE.out.winstats
    versions  = ch_versions

}