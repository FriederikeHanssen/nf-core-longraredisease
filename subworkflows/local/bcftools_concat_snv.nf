// Concatenate vcf files from clair3 and deepvariant if both are chosen


include { BCFTOOLS_CONCAT } from '../../modules/nf-core/bcftools/concat/main.nf'

workflow bcftools_concat_snv {
    take:
    input_vcfs  // channel: tuple(val(meta), path(vcfs), path(tbi))

    main:
    ch_versions = Channel.empty()
    BCFTOOLS_CONCAT(
        input_vcfs  // tuple(meta, vcfs, tbi)
    )

    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)
    emit:
    vcf = BCFTOOLS_CONCAT.out.vcf
    tbi = BCFTOOLS_CONCAT.out.tbi
    csi = BCFTOOLS_CONCAT.out.csi
    versions = ch_versions


}