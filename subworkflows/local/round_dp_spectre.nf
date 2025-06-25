include { ROUND_DP } from '../../modules/local/round_dp/main'

workflow round_dp_spectre_subworkflow {
    // Example input channel
    take:
    input_vcf

    main:
    
    ROUND_DP(input_vcf)

    emit:
    vcf = ROUND_DP.out.vcf
    versions = ROUND_DP.out.versions
}