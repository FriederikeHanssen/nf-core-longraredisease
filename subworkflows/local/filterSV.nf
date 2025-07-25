//
// Structural variant call filtering based on coverage and genomic regions
//

include { FILTERCOV_SV } from '../../modules/local/filterSV/main'

workflow filterSV_subworkflow {
    
    take:
    ch_sv_vcf           // channel: [meta, vcf] from SV callers (sniffles/cutesv/svim/survivor)
    ch_mosdepth_summary // channel: [meta, summary] from mosdepth_cnv_depth_subworkflow
    ch_bed_file         // channel: bed_file (target regions)
    chromosome_codes    // val: list of chromosome codes
    min_read_support    // val: minimum read support
    min_read_support_limit // val: minimum read support limit
    filter_pass        // val: boolean to filter PASS variants (optional, default: false)

    main:
    
    ch_versions = Channel.empty()
    
    // Filter SV calls using coverage information
    FILTERCOV_SV (
        ch_sv_vcf,
        ch_mosdepth_summary,
        ch_bed_file,
        chromosome_codes,
        min_read_support,
        min_read_support_limit,
        filter_pass
    )
    ch_versions = ch_versions.mix(FILTERCOV_SV.out.versions)
    
    emit:
    filtered_vcf = FILTERCOV_SV.out.filterbycov_vcf  
    versions     = ch_versions
}