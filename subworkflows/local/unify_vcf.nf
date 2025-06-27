//
// Unify multiple structural variant VCF files (SV, CNV, repeat)
//

include { UNIFY_VCFS } from '../../modules/local/unifyvcf_geneyx/main'

workflow unify_vcf_subworkflow {
    
    take:
    ch_sv_vcfs          // channel: [meta, [sv1.vcf, sv2.vcf, ...]] - Multiple SV VCF files
    ch_cnv_vcf          // channel: [meta, cnv.vcf] - Single CNV VCF file (optional)
    ch_repeat_vcf       // channel: [meta, repeat.vcf] - Single repeat VCF file (optional)
    
    main:
    
    ch_versions = Channel.empty()
    
    // Unify all VCF files
    UNIFY_VCFS(
        ch_sv_vcfs,
        ch_cnv_vcf,
        ch_repeat_vcf
    )
    ch_versions = ch_versions.mix(UNIFY_VCFS.out.versions)
    
    emit:
    unified_vcf = UNIFY_VCFS.out.unified_vcf    // channel: [meta, unified.vcf]
    versions    = ch_versions                   // channel: versions.yml
}