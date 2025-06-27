//
// Structural variant call filtering based on coverage and genomic regions
//

include { FILTERBYCOV_SV as FILTERBYCOV_SNIFFLES } from '../../modules/local/filterbycov_sv/main'

workflow filterbycov_sniffles_subworkflow {
    
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
    FILTERBYCOV_SNIFFLES(
        ch_sv_vcf,
        ch_mosdepth_summary,
        ch_bed_file,
        chromosome_codes,
        min_read_support,
        min_read_support_limit
    )
    ch_versions = ch_versions.mix(FILTERBYCOV_SNIFFLES.out.versions)
    
    emit:
    filtered_vcf = FILTERBYCOV_SNIFFLES.out.filterbycov_vcf  // ✅ Updated to match new emit name
    versions     = ch_versions
}