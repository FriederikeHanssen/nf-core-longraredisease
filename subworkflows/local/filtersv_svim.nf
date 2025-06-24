//
// Structural variant call filtering based on coverage and genomic regions
//

include { FILTER_SV_CALLS as FILTER_SV_CALL_SVIM } from '../../modules/local/filter_sv_calls'

workflow filterbycov_sv_svim {
    
    take:
    ch_sv_vcf           // channel: [meta, vcf] from SV callers (sniffles/cutesv/svim/survivor)
    ch_mosdepth_summary // channel: [meta, summary] from mosdepth_cnv_depth_subworkflow
    ch_bed_file         // channel: bed_file (target regions)
    chromosome_codes    // val: list of chromosome codes
    min_read_support    // val: minimum read support
    min_read_support_limit // val: minimum read support limit
    
    main:
    
    ch_versions = Channel.empty()
   
    
    
    // Filter SV calls using coverage information
    FILTER_SV_CALL_SVIM(
        ch_sv_vcf,
        ch_mosdepth_summary,
        ch_bed_file,
        chromosome_codes,
        min_read_support,
        min_read_support_limit
    )
    ch_versions = ch_versions.mix(FILTER_SV_CALL_SVIM.out.versions)

    emit:
    filtered_vcf = FILTER_SV_CALL_SVIM.out.filtered_vcf  // channel: [meta, filtered_vcf]
    versions     = ch_versions                        // channel: versions.yml
}