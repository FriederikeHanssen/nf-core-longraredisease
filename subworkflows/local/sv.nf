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
    input         // tuple(val(meta), path(bam), path(bai))
    fasta         // tuple(val(meta), path(fasta))
    tandem_file   // tuple(val(meta), path(bed))
    vcf_output    // val(true)
    snf_output    // val(true)
    primary_sv_caller // only relevant if consensuSV is false
    // filter sv inputs
    filter_sv            // boolean to filter SV calls based on coverage
    ch_mosdepth_summary // channel: [meta, summary] from mosdepth_cnv_depth_subworkflow
    ch_mosdepth_bed        // channel: bed_file (target regions)
    chromosome_codes    // val: list of chromosome codes
    min_read_support    // val: minimum read support
    min_read_support_limit // val: minimum read support limit
    filter_pass        // val: boolean to filter PASS variants (optional, default: false)
   
    main:
    // Run all SV callers
    SNIFFLES(input, fasta, tandem_file, vcf_output, snf_output)
    SVIM(input, fasta)
    CUTESV(input, fasta)

    // Sort and compress VCF files from SVIM
    BCFTOOLS_SORT_SVIM(SVIM.out.vcf)
    BGZIP_SVIM(BCFTOOLS_SORT_SVIM.out.vcf)  
    
    // Sort and compress VCF files from CUTESV
    BCFTOOLS_SORT_CUTESV(CUTESV.out.vcf)
    BGZIP_CUTESV(BCFTOOLS_SORT_CUTESV.out.vcf)

    if (filter_sv) {
        // Filter all SV callers
        ch_sniffles_vcf_gz_tbi = SNIFFLES.out.vcf
            .join(SNIFFLES.out.tbi, by: 0)
        
        FILTER_SV_SNIFFLES (
            ch_sniffles_vcf_gz_tbi,
            ch_mosdepth_summary,
            ch_mosdepth_bed,
            chromosome_codes,
            min_read_support,
            min_read_support_limit,
            filter_pass
        )

        ch_sniffles_vcf_gz = FILTER_SV_SNIFFLES.out.filterbycov_vcf
            .map { meta, vcf_gz, tbi -> 
                tuple(meta, vcf_gz) 
            }

        FILTER_SV_SVIM (
            BGZIP_SVIM.out.gz_tbi,
            ch_mosdepth_summary,
            ch_mosdepth_bed,
            chromosome_codes,
            min_read_support,
            min_read_support_limit,
            filter_pass
        )

        ch_svim_vcf_gz = FILTER_SV_SVIM.out.filterbycov_vcf
            .map { meta, vcf_gz, tbi -> 
                tuple(meta, vcf_gz) 
            }

        FILTER_SV_CUTESV (
            BGZIP_CUTESV.out.gz_tbi,
            ch_mosdepth_summary,
            ch_mosdepth_bed,
            chromosome_codes,
            min_read_support,
            min_read_support_limit,
            filter_pass
        )

        ch_cutesv_vcf_gz = FILTER_SV_CUTESV.out.filterbycov_vcf
            .map { meta, vcf_gz, tbi -> 
                tuple(meta, vcf_gz) 
            }
        
    } else {
        // No filtering - use original outputs
        ch_sniffles_vcf_gz = SNIFFLES.out.vcf
        ch_svim_vcf_gz = BGZIP_SVIM.out.gz_tbi
            .map { meta, vcf_gz, tbi -> 
                tuple(meta, vcf_gz) 
            }
        ch_cutesv_vcf_gz = BGZIP_CUTESV.out.gz_tbi
            .map { meta, vcf_gz, tbi -> 
                tuple(meta, vcf_gz) 
            }
    }

    // Select primary caller once (after filtering decision)
    if (primary_sv_caller == 'sniffles') {
        ch_vcf_gz = ch_sniffles_vcf_gz
    } else if (primary_sv_caller == 'cutesv') {
        ch_vcf_gz = ch_cutesv_vcf_gz
    } else if (primary_sv_caller == 'svim') {
        ch_vcf_gz = ch_svim_vcf_gz
    } else {
        log.warn "Unknown primary SV caller: ${primary_sv_caller}. Defaulting to Sniffles output."
        ch_vcf_gz = ch_sniffles_vcf_gz
    }

    emit:
    sniffles_vcf_gz      = ch_sniffles_vcf_gz
    sniffles_tbi         = SNIFFLES.out.tbi             
    sniffles_snf         = SNIFFLES.out.snf        
    svim_vcf             = BCFTOOLS_SORT_SVIM.out.vcf
    svim_vcf_gz          = ch_svim_vcf_gz  
    cutesv_vcf           = BCFTOOLS_SORT_CUTESV.out.vcf
    cutesv_vcf_gz        = ch_cutesv_vcf_gz   
    vcf_gz               = ch_vcf_gz
    sniffles_versions    = SNIFFLES.out.versions
    cutesv_versions      = CUTESV.out.versions
    svim_versions        = SVIM.out.versions   
}