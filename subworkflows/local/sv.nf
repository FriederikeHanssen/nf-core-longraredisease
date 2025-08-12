// Run sniffles SV calling
include { SNIFFLES                            } from '../../modules/nf-core/sniffles/main.nf'
// Run svim SV calling
include { SVIM                                } from '../../modules/local/svim/main.nf'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_SVIM } from '../../modules/nf-core/bcftools/sort/main.nf'
include { TABIX_BGZIPTABIX as BGZIP_SVIM      } from '../../modules/nf-core/tabix/bgziptabix/main.nf'
// Run cutesv SV calling
include { CUTESV                                } from '../../modules/nf-core/cutesv/main.nf'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_CUTESV } from '../../modules/nf-core/bcftools/sort/main.nf'
include { TABIX_BGZIPTABIX as BGZIP_CUTESV      } from '../../modules/nf-core/tabix/bgziptabix/main.nf'
include { FILTERCOV_SV as FILTER_SV_SNIFFLES  } from '../../modules/local/filterSV/main'
include { FILTERCOV_SV as FILTER_SV_SVIM      } from '../../modules/local/filterSV/main'
include { FILTERCOV_SV as FILTER_SV_CUTESV    } from '../../modules/local/filterSV/main'

workflow sv_subworkflow {
    
    take:
    input                    // tuple(val(meta), path(bam), path(bai))
    fasta                    // tuple(val(meta), path(fasta))
    tandem_file              // tuple(val(meta), path(bed))
    vcf_output               // val(true)
    snf_output               // val(true)
    primary_sv_caller        // val: primary caller name
    filter_sv                // val: boolean to filter SV calls
    ch_mosdepth_summary      // channel: [meta, summary]
    ch_mosdepth_bed          // channel: bed_file (target regions)
    chromosome_codes         // val: list of chromosome codes
    min_read_support         // val: minimum read support
    min_read_support_limit   // val: minimum read support limit
    filter_pass              // val: boolean to filter PASS variants
   
    main:
    ch_versions = Channel.empty()

    // Run all SV callers
    SNIFFLES(input, fasta, tandem_file, vcf_output, snf_output)
    SVIM(input, fasta)
    CUTESV(input, fasta)
    
    ch_versions = ch_versions.mix(SNIFFLES.out.versions)
    ch_versions = ch_versions.mix(SVIM.out.versions)
    ch_versions = ch_versions.mix(CUTESV.out.versions)

    // Sort and compress VCF files from SVIM and CUTESV
    BCFTOOLS_SORT_SVIM(SVIM.out.vcf)
    BGZIP_SVIM(BCFTOOLS_SORT_SVIM.out.vcf)  
    
    BCFTOOLS_SORT_CUTESV(CUTESV.out.vcf)
    BGZIP_CUTESV(BCFTOOLS_SORT_CUTESV.out.vcf)

    // Prepare channels for filtering or direct use
    ch_sniffles_input = SNIFFLES.out.vcf.join(SNIFFLES.out.tbi, by: 0)
    ch_svim_input = BGZIP_SVIM.out.gz_tbi
    ch_cutesv_input = BGZIP_CUTESV.out.gz_tbi

    // Apply filtering if requested
    if (filter_sv) {
        FILTER_SV_SNIFFLES(
            ch_sniffles_input,
            ch_mosdepth_summary,
            ch_mosdepth_bed,
            chromosome_codes,
            min_read_support,
            min_read_support_limit,
            filter_pass
        )

        FILTER_SV_SVIM(
            ch_svim_input,
            ch_mosdepth_summary,
            ch_mosdepth_bed,
            chromosome_codes,
            min_read_support,
            min_read_support_limit,
            filter_pass
        )

        FILTER_SV_CUTESV(
            ch_cutesv_input,
            ch_mosdepth_summary,
            ch_mosdepth_bed,
            chromosome_codes,
            min_read_support,
            min_read_support_limit,
            filter_pass
        )

        // Use filtered outputs
        ch_sniffles_vcf_tbi = FILTER_SV_SNIFFLES.out.filterbycov_vcf
        ch_svim_vcf_tbi = FILTER_SV_SVIM.out.filterbycov_vcf
        ch_cutesv_vcf_tbi = FILTER_SV_CUTESV.out.filterbycov_vcf
        
    } else {
        // Use original outputs
        ch_sniffles_vcf_tbi = ch_sniffles_input
        ch_svim_vcf_tbi = ch_svim_input
        ch_cutesv_vcf_tbi = ch_cutesv_input
    }

    // Extract VCF and TBI separately for individual outputs
    ch_sniffles_vcf_gz = ch_sniffles_vcf_tbi.map { meta, vcf_gz, tbi -> tuple(meta, vcf_gz) }
    ch_sniffles_tbi = ch_sniffles_vcf_tbi.map { meta, vcf_gz, tbi -> tuple(meta, tbi) }
    
    ch_svim_vcf_gz = ch_svim_vcf_tbi.map { meta, vcf_gz, tbi -> tuple(meta, vcf_gz) }
    ch_svim_tbi = ch_svim_vcf_tbi.map { meta, vcf_gz, tbi -> tuple(meta, tbi) }
    
    ch_cutesv_vcf_gz = ch_cutesv_vcf_tbi.map { meta, vcf_gz, tbi -> tuple(meta, vcf_gz) }
    ch_cutesv_tbi = ch_cutesv_vcf_tbi.map { meta, vcf_gz, tbi -> tuple(meta, tbi) }

    // Select primary caller outputs using proper Groovy if-else
    if (primary_sv_caller == 'sniffles') {
        ch_primary_vcf_gz = ch_sniffles_vcf_gz
        ch_primary_tbi = ch_sniffles_tbi
    } else if (primary_sv_caller == 'cutesv') {
        ch_primary_vcf_gz = ch_cutesv_vcf_gz
        ch_primary_tbi = ch_cutesv_tbi
    } else if (primary_sv_caller == 'svim') {
        ch_primary_vcf_gz = ch_svim_vcf_gz
        ch_primary_tbi = ch_svim_tbi
    } else {
        log.warn "Unknown primary SV caller: ${primary_sv_caller}. Defaulting to Sniffles."
        ch_primary_vcf_gz = ch_sniffles_vcf_gz
        ch_primary_tbi = ch_sniffles_tbi
    }

    emit:
    // Individual caller outputs
    sniffles_vcf_gz      = ch_sniffles_vcf_gz
    sniffles_tbi         = ch_sniffles_tbi
    sniffles_snf         = SNIFFLES.out.snf
    
    svim_vcf_gz          = ch_svim_vcf_gz
    svim_tbi             = ch_svim_tbi
    svim_vcf             = BCFTOOLS_SORT_SVIM.out.vcf
    
    cutesv_vcf_gz        = ch_cutesv_vcf_gz
    cutesv_tbi           = ch_cutesv_tbi
    cutesv_vcf           = BCFTOOLS_SORT_CUTESV.out.vcf
    
    // Primary caller outputs
    primary_vcf_gz       = ch_primary_vcf_gz
    primary_tbi          = ch_primary_tbi
    
    // Legacy output for backward compatibility
    vcf_gz               = ch_primary_vcf_gz
    
    // Version information
    versions             = ch_versions
}