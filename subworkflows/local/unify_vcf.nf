//
// Unify multiple structural variant VCF files (SV, CNV, repeat)
//

include { UNIFYVCF } from '../../modules/local/unifyvcf_geneyx/main'

workflow unify_vcf_subworkflow {
    
    take:
    ch_sv_vcfs          // channel: [meta, [sv1.vcf, sv2.vcf, ...]] - Multiple SV VCF files
    ch_cnv_vcf          // channel: [meta, cnv.vcf] - Single CNV VCF file (optional)
    ch_repeat_vcf       // channel: [meta, repeat.vcf] - Single repeat VCF file (optional)
    modify_repeats      // Boolean: whether to modify repeat calls (true for STRaglr)

    main:
    
    ch_versions = Channel.empty()
    
    // Unify all VCF files
    UNIFYVCF(
        ch_sv_vcfs,
        ch_cnv_vcf,
        ch_repeat_vcf,
        modify_repeats
    )
    ch_versions = ch_versions.mix(UNIFYVCF.out.versions)
    
    emit:
    unified_vcf = UNIFYVCF.out.unified_vcf    // channel: [meta, unified.vcf]
    unified_tbi = UNIFYVCF.out.unified_tbi    // channel: [meta, unified.vcf.tbi]
    versions    = ch_versions                   // channel: versions.yml
}