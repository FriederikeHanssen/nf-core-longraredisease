#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NANORAREDX: Comprehensive Nanopore Rare Disease Analysis Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    This pipeline performs comprehensive genomic analysis of nanopore sequencing data
    for rare disease applications, including:
    - Structural variant (SV) calling with multiple callers
    - Single nucleotide variant (SNV) calling  
    - Copy number variant (CNV) detection
    - Short tandem repeat (STR) analysis
    - Phasing analysis
    - Optional methylation calling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Data preprocessing subworkflows
include { BAM_STATS_SAMTOOLS                 } from '../subworkflows/nf-core/bam_stats_samtools/main.nf'
include { samtools_index_subworkflow         } from '../subworkflows/local/samtools_index.nf'
include { generate_fai_subworkflow           } from '../subworkflows/local/generate_fai.nf'
include { samtools_merge_bam_subworkflow     } from '../subworkflows/local/samtools_merge_bam.nf'
include { samtools_bam_to_fastq_subworkflow  } from '../subworkflows/local/samtools_bam_to_fastq.nf'
include { minimap2_align_subworkflow         } from '../subworkflows/local/minimap2_align.nf'
include { nanoplot_subworkflow               } from '../subworkflows/local/nanoplot.nf'

// Methylation calling 
include { modkit_methyl_subworkflow } from '../subworkflows/local/modkit_methyl.nf'
include { modkit_cpg_subworkflow    } from '../subworkflows/local/modkit_cpg.nf'

// Coverage analysis subworkflows
include { mosdepth_subworkflow              } from '../subworkflows/local/mosdepth.nf'
include { mosdepth_cov_analysis_subworkflow } from '../subworkflows/local/mosdepth_cov_analysis.nf'

// Structural variant calling subworkflows
include { sniffles_sv_subworkflow } from '../subworkflows/local/sniffles_sv.nf'
include { cutesv_sv_subworkflow   } from '../subworkflows/local/cutesv_sv.nf'
include { svim_sv_subworkflow     } from '../subworkflows/local/svim_sv.nf'

// SV filtering subworkflows - coverage-based filtering
include { low_coverage_sv_filter_subworkflow } from '../subworkflows/local/low_coverage_sv_filter.nf'

// SV merging and intersection filtering subworkflows
include { multi_caller_sv_filter_subworkflow } from '../subworkflows/local/multi_caller_sv_filter.nf'

// SNV calling and processing subworkflows
include { clair3_snv_subworkflow          } from '../subworkflows/local/clair3_snv.nf'
include { deepvariant_snv_subworkflow     } from '../subworkflows/local/deepvariant_snv.nf'
include { bcftools_concat_snv_subworkflow } from '../subworkflows/local/bcftools_concat_snv.nf'

// Filter SNV VCF by PASS status
include { bcftools_filter_clair3_subworkflow } from '../subworkflows/local/bcftools_filter_clair3.nf'
include { bcftools_filter_deepvariant_subworkflow } from '../subworkflows/local/bcftools_filter_deepvariant.nf'

// Phasing subworkflow
include { longphase_subworkflow } from '../subworkflows/local/longphase.nf'

// CNV calling subworkflows
include { spectre_cnv_subworkflow      } from '../subworkflows/local/spectre_cnv.nf'
include { round_dp_spectre_subworkflow } from '../subworkflows/local/round_dp_spectre.nf'
include { qdnaseq_cnv_subworkflow      } from '../subworkflows/local/qdnaseq_cnv.nf'

// STR analysis subworkflow
include { straglr_str_subworkflow } from '../subworkflows/local/straglr_str.nf'

// VCF processing subworkflows
include { unify_vcf_subworkflow } from '../subworkflows/local/unify_vcf.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow nanoraredx {

    // Parameter validation
    if (params.phase && !params.snv) {
        error "Phasing requires SNV calling to be enabled (--snv true)"
    }

    if (params.phase_with_sv && !params.sv) {
        error "Phasing with SVs requires SV calling to be enabled (--sv true)"
    }

    if (params.exclude_low_coverage_SV && !params.sv) {
        error "SV coverage filtering requires SV calling to be enabled (--sv true)"
    }

    if (params.unify_SV_CNV_STR_vcfs && !(params.sv && params.cnv && params.str)) {
        error "VCF unification requires SV, CNV, and STR to be enabled (--sv true, --cnv true, --str true)"
    }

    if (params.cnv && !params.run_qdnaseq && !params.use_test_data && !params.snv) {
        error "Spectre CNV calling requires SNV calling to be enabled (--snv true) unless using test data"
    }

    if (params.multi_caller_sv_filtering && !params.sv) {
        error "Multi-caller SV filtering requires SV calling to be enabled (--sv true)"
    }

    // Input validation based on alignment parameter
    if (params.align) {
        // When align=true, expect unaligned BAMs in bam_dir
        if (!params.bam_dir) {
            error "When --align is true, --bam_dir must be provided with unaligned BAM files"
        }
        if (!params.sample_name) {
            error "When --align is true, --sample_name must be provided"
        }
    } else {
        // When align=false, expect aligned BAM file
        if (!params.aligned_bam) {
            error "When --align is false, --aligned_bam must be provided with the path to aligned BAM file"
        }
        // Extract sample name from aligned BAM filename if not provided
        if (!params.sample_name) {
            def aligned_bam_file = file(params.aligned_bam)
            // Handle naming pattern like "A.sorted.bam" -> extract "A"
            def extracted_sample_name = aligned_bam_file.name
                .replaceAll(/\.sorted\.bam$/, '')  // Remove .sorted.bam
                .replaceAll(/\.bam$/, '')          // Remove .bam (fallback)
            
            if (extracted_sample_name.isEmpty()) {
                error "Could not extract sample name from aligned BAM filename: ${aligned_bam_file.name}. Please provide --sample_name manually."
            }
            
            params.sample_name = extracted_sample_name
            log.info "Sample name extracted from aligned BAM: ${params.sample_name}"
        }
    }

/*
=======================================================================================
                                REFERENCE FILES SETUP
=======================================================================================
*/
    
    ch_fasta = Channel
    .fromPath(params.fasta_file, checkIfExists: true)
    .map { fasta -> tuple([id: "ref"], fasta) }

    // Generate FAI index
    generate_fai_subworkflow(ch_fasta, true)
    ch_fai = generate_fai_subworkflow.out.fai

    // Create combined FASTA+FAI channel by joining
    ch_fasta_fai = ch_fasta
    .join(ch_fai, by: 0)
    .map { meta, fasta, fai -> tuple(meta, fasta, fai) }
    

    // Tandem repeat file for Sniffles (only if SV calling is enabled)
    if (params.sv) {
        ch_trf = Channel
            .fromPath(params.sniffles_tandem_file, checkIfExists: true)
            .map { bed -> tuple([id: "trf"], bed) }
    } else {
        ch_trf = Channel.empty()
    }

/*
=======================================================================================
                               DATA PREPROCESSING PIPELINE
=======================================================================================
*/
    
    if (params.align) {
        /*
        ================================================================================
                            ALIGNMENT WORKFLOW (UNALIGNED INPUT)
        ================================================================================
        */
        
        log.info "Sample name: ${params.sample_name}"
        
        // Collect unaligned BAM files and use parameter-defined sample name
        ch_bam_files = Channel
            .fromPath("${params.bam_dir}/*.bam")
            .map { bam ->
                // Use the sample_name parameter for unaligned BAMs
                return [ [id: params.sample_name], bam ]
            }
            .groupTuple()

        // Merge multiple unaligned BAM files per sample into single BAM
        samtools_merge_bam_subworkflow(
            ch_bam_files, 
            [[:], []], 
            [[:], []]
        )

        // Convert merged BAMs to FASTQ format while preserving methylation tags
        samtools_bam_to_fastq_subworkflow(
            samtools_merge_bam_subworkflow.out.bam
        )

        // Align FASTQ reads to reference genome using minimap2
        minimap2_align_subworkflow(
            ch_fasta,
            samtools_bam_to_fastq_subworkflow.out.other
        )

        // Set final aligned BAM channels from minimap2 output
        ch_final_sorted_bam = minimap2_align_subworkflow.out.ch_sorted_bam
        ch_final_sorted_bai = minimap2_align_subworkflow.out.ch_sorted_bai

        // Prepare input for nanoplot from FASTQ
        ch_fastq_nanoplot = samtools_bam_to_fastq_subworkflow.out.other
            .map { meta, fastq_file ->
                tuple(meta, fastq_file) 
            }

    } else {
        /*
        ================================================================================
                            ALIGNED INPUT WORKFLOW (ALIGNED BAM INPUT)
        ================================================================================
        */
        
        // Use pre-aligned BAM file
        ch_aligned_bam = Channel
            .fromPath(params.aligned_bam, checkIfExists: true)
            .map { bam ->
                def meta = [id: params.sample_name]
                return [meta, bam]
            }

        // Check if BAI file already exists before generating
        def bai_file = file("${params.aligned_bam}.bai")
        
        if (bai_file.exists()) {
            // Create channel with existing BAI file
            ch_aligned_bai = ch_aligned_bam.map { meta, bam ->
                return [meta, file("${bam}.bai")]
            }
        } else {
            // Generate BAI index through proper subworkflow
            samtools_index_subworkflow(ch_aligned_bam)
            ch_aligned_bai = samtools_index_subworkflow.out.bai
        }

        // Set final aligned BAM channels from input
        ch_final_sorted_bam = ch_aligned_bam
        ch_final_sorted_bai = ch_aligned_bai

        // For nanoplot, we'll skip it (no FASTQ available)
        ch_fastq_nanoplot = Channel.empty()
    }

/*
=======================================================================================
                                COVERAGE ANALYSIS
=======================================================================================
*/
    
    // Prepare input channel with BAM, BAI, and optional BED file for coverage analysis
    ch_input_bam_bai_bed = ch_final_sorted_bam
        .join(ch_final_sorted_bai, by: 0)
        .map { meta, bam, bai ->
            def bed = params.bed_file ? file(params.bed_file) : []
            tuple(meta, bam, bai, bed)
        }

    // Prepare simplified BAM input channel for variant calling and methylation calling
    ch_input_bam = ch_final_sorted_bam
        .join(ch_final_sorted_bai, by: 0)
        .map { meta, bam, bai -> tuple(meta, bam, bai) }

    if (params.generate_bam_stats) {
        // Run BAM statistics using samtools
        BAM_STATS_SAMTOOLS(
            ch_input_bam,
            ch_fasta
        )
    }
    
    // Run nanoplot (only if we have FASTQ data from alignment workflow)
    if (params.qc && params.align) {
        nanoplot_subworkflow(
            ch_fastq_nanoplot
        )
    } else if (params.qc && !params.align) {
        log.warn "Skipping NanoPlot: QC enabled but no FASTQ data available (using pre-aligned BAM)"
    }
    
    // Run mosdepth when needed
    if (params.sv || params.cnv || params.generate_coverage) {
        mosdepth_subworkflow(
            ch_input_bam_bai_bed,
            [[:], []]
        )
    }

    // Coverage analysis only needed for SV filtering
    if (params.sv && params.exclude_low_coverage_SV) {
        mosdepth_cov_analysis_subworkflow(
            ch_input_bam_bai_bed,
            [[:], []]
        )
    }

    // Methylation calling with modkit (if enabled)
    ch_empty_bed = Channel.value([[:], []])

    if (params.methyl) {
        modkit_methyl_subworkflow(
            ch_input_bam,
            ch_fasta_fai,
            ch_empty_bed  // Empty channel for optional BED file
        )

        modkit_cpg_subworkflow(
            ch_input_bam,
            ch_fasta_fai,
            ch_empty_bed 
        )
    }

   
/*
=======================================================================================
                        STRUCTURAL VARIANT CALLING WORKFLOW
=======================================================================================
*/

if (params.sv) {

    /*
    ================================================================================
                            PARALLEL SV CALLER EXECUTION
    ================================================================================
    */
    
    // Sniffles: Population-scale SV calling (PRIMARY CALLER)
    sniffles_sv_subworkflow(
        ch_input_bam,
        ch_fasta,
        ch_trf,
        params.vcf_output,
        params.snf_output
    )

    // CuteSV: Long-read SV detection with high sensitivity
    cutesv_sv_subworkflow(
        ch_input_bam,
        ch_fasta
    )

    // SVIM: Structural variant identification using alignment information
    svim_sv_subworkflow(
        ch_input_bam,
        ch_fasta
    )

    // Collect outputs from all SV callers
    ch_sniffles_vcf = sniffles_sv_subworkflow.out.vcf
    ch_svim_vcf     = svim_sv_subworkflow.out.vcf
    ch_cutesv_vcf   = cutesv_sv_subworkflow.out.vcf
    
    // Initialize primary SV VCF with Sniffles output (our gold standard)
    ch_primary_sv_vcf = ch_sniffles_vcf

    /*
    ================================================================================
                        MULTI-CALLER FILTERING AND CONSENSUS
    ================================================================================
    */
    
    if (params.multi_caller_sv_filtering) {
        // Create chromosome sizes channel
        ch_chrom_sizes = params.chrom_sizes ? 
            Channel.value([[:], file(params.chrom_sizes, checkIfExists: true)]) : 
            Channel.value([[:], []])

        // Prepare VCFs for SURVIVOR merging
        ch_vcfs_for_merging = ch_sniffles_vcf
            .map { meta, vcf -> 
                [params.sample_name ?: meta.id, vcf, 'sniffles'] 
            }
            .mix(
                ch_cutesv_vcf
                    .filter { meta, vcf -> vcf && vcf.exists() }
                    .map { meta, vcf -> 
                        [params.sample_name ?: meta.id, vcf, 'cutesv'] 
                    }
            )
            .mix(
                ch_svim_vcf
                    .filter { meta, vcf -> vcf && vcf.exists() }
                    .map { meta, vcf -> 
                        [params.sample_name ?: meta.id, vcf, 'svim'] 
                    }
            )
            .groupTuple(by: 0)
            .map { sample_id, vcfs, callers ->
                def meta = [id: sample_id, callers: callers]
                [meta, vcfs]
            }

        // Run multi-caller filtering
        multi_caller_sv_filter_subworkflow(
            ch_vcfs_for_merging,
            ch_primary_sv_vcf,                                // Use consistent channel name
            ch_chrom_sizes,
            params.max_distance_breakpoints ?: 1000,
            params.min_supporting_callers ?: 2,
            params.account_for_type ?: true,
            params.account_for_sv_strands ?: true,
            params.estimate_distanced_by_sv_size ?: false,
            params.min_sv_size ?: 50
        )

        // Update primary VCF to filtered result
        ch_primary_sv_vcf = multi_caller_sv_filter_subworkflow.out.filtered_vcf_sniffles
    }

    /*
    ================================================================================
                            COVERAGE-BASED FILTERING
    ================================================================================
    */
    
    if (params.exclude_low_coverage_SV) {
        
        // Apply coverage-based filtering
        low_coverage_sv_filter_subworkflow(
            ch_primary_sv_vcf,                                // Use consistent channel name
            mosdepth_subworkflow.out.summary_txt,
            mosdepth_cov_analysis_subworkflow.out.quantized_bed,
            params.chromosome_codes ?: 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY',
            params.min_read_support ?: 'auto',
            params.min_read_support_limit ?: 3,
            params.filter_pass ?: false
        )
        
        // Update primary VCF to coverage-filtered result
        ch_primary_sv_vcf = low_coverage_sv_filter_subworkflow.out.filtered_vcf
    }
 
    } 
    
else {
    
    // Create empty channels for consistent workflow structure
    ch_sniffles_vcf   = Channel.empty()
    ch_svim_vcf       = Channel.empty()
    ch_cutesv_vcf     = Channel.empty()
    ch_primary_sv_vcf = Channel.empty()
}
   
/*
================================================================================
                        SINGLE NUCLEOTIDE VARIANT CALLING
================================================================================
*/

if (params.snv) {

    
    // Prepare input for Clair3 SNV calling
    ch_input_bam_clair3 = ch_input_bam.map { meta, bam, bai ->
        def id = bam.getBaseName().replaceAll(/\.bam$/, '')
        def updated_meta = meta + [id: id]
        tuple(
            updated_meta,
            bam,
            bai,
            params.clair3_model,
            [],
            params.clair3_platform
        )
    }
    
    // Run Clair3 SNV calling
    clair3_snv_subworkflow(
        ch_input_bam_clair3,
        ch_fasta,
        ch_fai
    )

    ch_clair3_vcf = clair3_snv_subworkflow.out.vcf
    ch_clair3_tbi = clair3_snv_subworkflow.out.tbi

    if (params.filter_pass) {
        // Prepare input channels for BCFTOOLS_VIEW
        ch_clair3_vcf_tbi = ch_clair3_vcf.join(ch_clair3_tbi, by: 0)
        
        // Create empty channels for optional inputs
        ch_regions = Channel.value([])
        ch_targets = Channel.value([])  
        ch_samples = Channel.value([])
        
        // Call BCFTOOLS_VIEW with correct input structure
        bcftools_filter_clair3_subworkflow(
            ch_clair3_vcf_tbi,   // tuple val(meta), path(vcf), path(index)
            ch_regions,          // path(regions) - empty
            ch_targets,          // path(targets) - empty  
            ch_samples           // path(samples) - empty
        )
        
        // Use filtered results
        ch_clair3_vcf = bcftools_filter_clair3_subworkflow.out.vcf
        ch_clair3_tbi = bcftools_filter_clair3_subworkflow.out.tbi
    } 
    } else {
    // No SNV calling - create empty channels
    ch_clair3_vcf = Channel.empty()
    ch_clair3_tbi = Channel.empty()
    }

// Run DeepVariant SNV calling if requested

if (params.run_deepvariant) {

    
    deepvariant_snv_subworkflow(
        ch_input_bam_bai_bed,
        ch_fasta,
        ch_fai,
        [[:], []],
        [[:], []]
    )
    
    ch_dv_vcf = deepvariant_snv_subworkflow.out.vcf
    ch_dv_tbi = deepvariant_snv_subworkflow.out.tbi

    if (params.filter_pass) {
        // Prepare input channels for BCFTOOLS_VIEW
        ch_dv_vcf_tbi = ch_dv_vcf.join(ch_dv_tbi, by: 0)
        
        // Create empty channels for optional inputs
        ch_regions = Channel.value([])
        ch_targets = Channel.value([])  
        ch_samples = Channel.value([])
        
        // Call BCFTOOLS_VIEW with correct input structure
        bcftools_filter_deepvariant_subworkflow(
            ch_dv_vcf_tbi,       // tuple val(meta), path(vcf), path(index)
            ch_regions,          // path(regions) - empty
            ch_targets,          // path(targets) - empty  
            ch_samples           // path(samples) - empty
        )
        
    } 
    }

    // Unify SNV VCFs from Clair3 and DeepVariant if requested
if (params.unify_SNV && params.run_deepvariant && params.snv) {

    
    // Combine and concatenate VCF files from both callers
    // The bcftools_concat_snv_subworkflow expects: tuple(meta, [vcfs], [tbis])
    combined_vcfs = ch_clair3_vcf
        .join(ch_clair3_tbi, by: 0)
        .join(ch_dv_vcf.join(ch_dv_tbi, by: 0), by: 0)
        .map { meta, clair3_vcf, clair3_tbi, dv_vcf, dv_tbi ->
            [
                meta,
                [clair3_vcf, dv_vcf],    // List of VCF files
                [clair3_tbi, dv_tbi]     // List of corresponding TBI files
            ]
        }

    // Run concatenation (output not captured since not needed elsewhere)
    bcftools_concat_snv_subworkflow(combined_vcfs)
    
    } 
  

/*
=======================================================================================
                                 PHASING ANALYSIS
=======================================================================================
*/

// Run phasing with LongPhase if enabled
if (params.phase && params.snv) {

    
    if (params.sv && params.phase_with_sv) {
        // Phasing with both SNVs and SVs
        ch_longphase_input = ch_input_bam
            .join(ch_clair3_vcf, by: 0)
            .join(ch_primary_sv_vcf, by: 0)
            .map { meta, bam, bai, snv_vcf, sv_vcf -> 
                tuple(meta, bam, bai, snv_vcf, sv_vcf, []) 
            }
            
            } 
    else {
        // Phasing with SNVs only
        ch_longphase_input = ch_input_bam
            .join(ch_clair3_vcf, by: 0)
            .map { meta, bam, bai, snv_vcf -> 
            tuple(meta, bam, bai, snv_vcf, [], []) 
            }   
            }

    longphase_subworkflow(
        ch_longphase_input,
        ch_fasta,
        ch_fai
    )

    }

/*
=======================================================================================
                        COPY NUMBER VARIANT CALLING
=======================================================================================
*/

if (params.cnv) {

     // Spectre CNV calling - requires SNV data unless using test data
        
        if (params.use_test_data) {

            // Test mode with hardcoded parameters as can not be run on subset of data
            ch_spectre_test_reference = Channel
                .fromPath(params.spectre_test_clair3_vcf, checkIfExists: true)
                .map { vcf_file -> 
                    return [id: params.sample_name]  // Use sample_name parameter
                }
                .combine(Channel.fromPath(params.spectre_test_fasta_file, checkIfExists: true))
                .map { meta, fasta -> tuple(meta, fasta) }
            
            spectre_cnv_subworkflow(
                params.spectre_test_mosdepth,
                ch_spectre_test_reference,
                params.spectre_test_clair3_vcf,
                params.spectre_metadata,
                params.spectre_blacklist
            )
            ch_spectre_vcf = spectre_cnv_subworkflow.out.vcf
            } 
            
        else {
            
            if (!params.snv) {
                error "Spectre CNV calling requires SNV calling to be enabled (--snv true) or use QDNAseq (--run_qdnaseq true)"
            }
            
            // Extract just the BED file from mosdepth output
            ch_spectre_mosdepth_bed = mosdepth_subworkflow.out.regions_bed.map { meta, bed -> bed }
            
            // Prepare reference channel for Spectre CNV calling
            ch_spectre_reference = ch_clair3_vcf
            .map { meta, vcf_file -> 
            [id: params.sample_name]  // Create new meta with sample_name
            }
            .combine(Channel.fromPath(params.fasta_file, checkIfExists: true))
            .map { meta, fasta -> tuple(meta, fasta) }
            
            // Extract just the VCF file from the ch_final_snv_vcf channel
            ch_spectre_clair3_vcf = ch_clair3_vcf.map { meta, vcf -> vcf }
    
            spectre_cnv_subworkflow(
                ch_spectre_mosdepth_bed,
                ch_spectre_reference, 
                ch_spectre_clair3_vcf,              
                params.spectre_metadata,
                params.spectre_blacklist
            )

            ch_spectre_vcf = spectre_cnv_subworkflow.out.vcf
        
         }

        // Round decimal places in Spectre output
        round_dp_spectre_subworkflow(ch_spectre_vcf)
        ch_cnv_vcf = round_dp_spectre_subworkflow.out.vcf
    }
    else {
   
    ch_cnv_vcf = Channel.empty()
    
    }

if (params.run_qdnaseq) {
        // QDNAseq doesn't need mosdepth or SNV data
        qdnaseq_cnv_subworkflow(
            ch_input_bam,
            params.genome_build,
            params.qdnaseq_bin_size,
            params.method,
            params.cutoff,
            params.cutoff_del,
            params.cutoff_loss,
            params.cutoff_gain,
            params.cellularity
        )
        }
/*
================================================================================
                        SHORT TANDEM REPEAT ANALYSIS
================================================================================
*/

if (params.str) {
    straglr_str_subworkflow(
        ch_input_bam,
        ch_fasta,
        params.str_bed_file
    )
    ch_str_vcf = straglr_str_subworkflow.out.vcf
    
    } 
    
    else {
    ch_str_vcf = Channel.empty()
    }

/*
================================================================================
                             VCF UNIFICATION 
================================================================================
*/

if (params.unify_SV_CNV_STR_vcfs) {
    if (params.sv && params.cnv && params.str) {
        
        // Prepare individual VCF channels with consistent metadata
        ch_sv_for_unify = ch_primary_sv_vcf.map { meta, vcf -> [[id: params.sample_name], vcf] }
        
        ch_cnv_for_unify = ch_cnv_vcf.map { meta, vcf -> [[id: params.sample_name], vcf] }
        
        ch_str_for_unify = ch_str_vcf.map { meta, vcf -> [[id: params.sample_name], vcf] }
        
        // Run VCF unification subworkflow
        unify_vcf_subworkflow(
            ch_sv_for_unify,                    // SV VCFs
            ch_cnv_for_unify,                   // CNV VCFs
            ch_str_for_unify,                   // STR VCFs
            params.modify_str_calls             // Modify repeat calls flag
        )
    } else {
        error "VCF unification requires all three variant calling types to be enabled: --sv true --cnv true --str true"
    }
    } 
    /*
    ================================================================================
                            FUTURE ENHANCEMENTS (COMMENTED OUT)
    ================================================================================
    */
    
    // The following sections are prepared for future implementation:
    
    // Methylation calling with MEOW
    // SV annotation with Savanna
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


