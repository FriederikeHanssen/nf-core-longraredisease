// Concatenate vcf files from clair3 and deepvariant if both are chosen


include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_SNV } from '../../modules/nf-core/bcftools/concat/main.nf'

workflow bcftools_concat_snv_subworkflow {
    take:
    input_vcfs  // channel: tuple(val(meta), path(vcfs), path(tbi))

    main:
    ch_versions = Channel.empty()
    BCFTOOLS_CONCAT_SNV(
        input_vcfs  // tuple(meta, vcfs, tbi)
    )

    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT_SNV.out.versions)
    emit:
    vcf = BCFTOOLS_CONCAT_SNV.out.vcf
    tbi = BCFTOOLS_CONCAT_SNV.out.tbi
    csi = BCFTOOLS_CONCAT_SNV.out.csi
    versions = ch_versions


}