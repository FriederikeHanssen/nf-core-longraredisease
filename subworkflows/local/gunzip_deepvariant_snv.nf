// Gunzip if deepvariant is run

include { GUNZIP as GUNZIP_DEEPVARIANT } from '../../modules/nf-core/gunzip/main.nf'

workflow gunczip_deepvariant_snv_subworkflow { 
    take:
    input_vcf  // channel: tuple(val(meta), path(vcf), path(tbi))

    main:
    ch_versions = Channel.empty()

    // Run GUNZIP on DeepVariant VCF files
    GUNZIP_DEEPVARIANT(input_vcf)

    ch_versions = ch_versions.mix(GUNZIP_DEEPVARIANT.out.versions)

    emit:
    vcf         = GUNZIP_DEEPVARIANT.out.gunzip  // Uncompressed VCF files
    versions    = ch_versions                     // Version information
}