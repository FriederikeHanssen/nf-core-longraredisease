include { BCFTOOLS_VIEW as BCFTOOLS_FILTER_DEEPVARIANT } from '../../modules/nf-core/bcftools/view/main.nf'

workflow bcftools_filter_deepvariant_subworkflow {

take:
    vcf
    regions
    targets
    samples

main:

    BCFTOOLS_FILTER_DEEPVARIANT(
        vcf,
        regions,
        targets,
        samples
    )


emit:

vcf  =  BCFTOOLS_FILTER_DEEPVARIANT.out.vcf
tbi  =  BCFTOOLS_FILTER_DEEPVARIANT.out.tbi
}