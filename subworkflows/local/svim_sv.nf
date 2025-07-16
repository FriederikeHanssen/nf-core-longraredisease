// SVIM structural variant caller - used a support caller

include { SVIM } from '../../modules/local/svim/main.nf'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_SVIM } from '../../modules/nf-core/bcftools/sort/main.nf'

workflow svim_sv_subworkflow {

    take:
    input         // tuple(val(meta), path(bam), path(bai))
    fasta         // tuple(val(meta), path(fasta))

    main:
    // 1. Call SVs with SVIM
    SVIM(input, fasta)

    // 2. Sort the VCF with BCFTOOLS_SORT
    BCFTOOLS_SORT_SVIM(SVIM.out.vcf)


    emit:
    vcf         = BCFTOOLS_SORT_SVIM.out.vcf
    versions    = SVIM.out.versions
}
