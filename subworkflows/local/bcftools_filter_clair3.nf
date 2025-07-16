include { BCFTOOLS_VIEW as BCFTOOLS_FILTER_CLAIR3 } from '../../modules/nf-core/bcftools/view/main.nf'

workflow bcftools_filter_clair3_subworkflow {

take:
    vcf
    regions
    targets
    samples

main:

    BCFTOOLS_FILTER_CLAIR3(
        vcf,
        regions,
        targets,
        samples
    )


emit:

vcf  =  BCFTOOLS_FILTER_CLAIR3.out.vcf
tbi  =  BCFTOOLS_FILTER_CLAIR3.out.tbi
}