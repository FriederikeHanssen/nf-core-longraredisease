// Concatenate vcf files from clair3 and deepvariant if both are chosen
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_SNV } from '../../modules/nf-core/bcftools/concat/main.nf'
// add bcftools isec

workflow merge_snv_subworkflow {
    
    take:
    consensuSNV_vcf  // channel: tuple(val(meta), path(vcfs), path(tbi))
   
    

    main:
   
    BCFTOOLS_CONCAT_SNV(
        consensuSNV_vcf   // tuple(meta, vcfs, tbi)
    )


    ch_versions = BCFTOOLS_CONCAT_SNV.out.versions
    
    emit:
    vcf = BCFTOOLS_CONCAT_SNV.out.vcf
    tbi = BCFTOOLS_CONCAT_SNV.out.tbi
    csi = BCFTOOLS_CONCAT_SNV.out.csi
    versions = ch_versions


}
