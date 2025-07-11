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

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Data preprocessing subworkflows
include { bam_stats_subworkflow              } from '../subworkflows/local/samtools_bam_stats.nf'
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
include { filterbycov_svim_subworkflow     } from '../subworkflows/local/filterbycov_svim.nf'
include { filterbycov_sniffles_subworkflow } from '../subworkflows/local/filterbycov_sniffles.nf'
include { filterbycov_cutesv_subworkflow   } from '../subworkflows/local/filterbycov_cutesv.nf'

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

    if (params.filter_sv_coverage && !params.sv) {
        error "SV coverage filtering requires SV calling to be enabled (--sv true)"
    }

    if (params.unify_vcfs && !(params.sv && params.cnv && params.str)) {
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
    ================================================================================
                                REFERENCE FILES SETUP
    ================================================================================
    */
    
    // Reference genome FASTA file
    ch_fasta = Channel
        .fromPath(params.fasta_file, checkIfExists: true)
        .map { fasta -> tuple([id: "ref"], fasta) }

    // Reference genome index file (generated from samtools faidx)
    def fai_file = file("${params.fasta_file}.fai")
    
    if (fai_file.exists()) {
        log.info "Using existing FAI file: ${fai_file}"
        ch_fai = Channel
            .fromPath("${params.fasta_file}.fai", checkIfExists: true)
            .map { fai -> tuple([id: "ref_fai"], fai) }
        
        ch_fasta_fai = Channel
            .fromPath(params.fasta_file, checkIfExists: true)
            .map { fasta -> 
                def fai = file("${fasta}.fai")
                tuple([id: "ref"], fasta, fai)
            }
    } else {
        log.info "FAI file not found. Generating index for: ${params.fasta_file}"
        
        // Generate FAI file using samtools faidx
        generate_fai_subworkflow(ch_fasta, true)
        
        ch_fai = generate_fai_subworkflow.out.ch_fai
        ch_fasta_fai = ch_fasta
            .join(generate_fai_subworkflow.out.ch_fai, by: 0)
            .map { meta, fasta, fai -> tuple(meta, fasta, fai) }
    }

    // Tandem repeat file for Sniffles (only if SV calling is enabled)
    if (params.sv) {
        ch_trf = Channel
            .fromPath(params.sniffles_tandem_file, checkIfExists: true)
            .map { bed -> tuple([id: "trf"], bed) }
    } else {
        ch_trf = Channel.empty()
    }

    /*
    ================================================================================
                            DATA PREPROCESSING PIPELINE
    ================================================================================
    */
    
    if (params.align) {
        /*
        ================================================================================
                            ALIGNMENT WORKFLOW (UNALIGNED INPUT)
        ================================================================================
        */
        
        log.info "Running alignment workflow with unaligned BAMs from: ${params.bam_dir}"
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
        
        log.info "Using pre-aligned BAM: ${params.aligned_bam}"
        log.info "Sample name: ${params.sample_name}"
        log.info "Skipping alignment steps (merge, fastq conversion, minimap2)"
        
        // Use pre-aligned BAM file
        ch_aligned_bam = Channel
            .fromPath(params.aligned_bam, checkIfExists: true)
            .map { bam ->
                def meta = [id: params.sample_name]
                return [meta, bam]
            }

        // Check if BAI file exists and handle accordingly
        def bai_file = file("${params.aligned_bam}.bai")
        if (bai_file.exists()) {
            log.info "Using existing BAI file: ${bai_file}"
            ch_aligned_bai = Channel
                .fromPath("${params.aligned_bam}.bai", checkIfExists: true)
                .map { bai ->
                    def meta = [id: params.sample_name]
                    return [meta, bai]
                }
        } else {
            log.info "BAI file not found. Generating index for: ${params.aligned_bam}"
            
            // Generate BAI file using samtools index
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
    ================================================================================
                                COVERAGE ANALYSIS
    ================================================================================
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
        log.info "Generating BAM statistics using samtools"
        // Run BAM statistics using samtools
        bam_stats_subworkflow(
            ch_input_bam,
            ch_fasta
        )
    }
    
    // Run nanoplot (only if we have FASTQ data from alignment workflow)
    if (params.qc && params.align) {
        log.info "Running NanoPlot for quality control"
        nanoplot_subworkflow(
            ch_fastq_nanoplot
        )
    } else if (params.qc && !params.align) {
        log.info "Skipping NanoPlot: QC enabled but no FASTQ data available (using pre-aligned BAM)"
    }
    
    // Run mosdepth when needed
    if (params.sv || params.cnv || params.generate_coverage) {
        log.info "Running mosdepth for coverage analysis"
        mosdepth_subworkflow(
            ch_input_bam_bai_bed,
            [[:], []]
        )
    }

    // Coverage analysis only needed for SV filtering
    if (params.sv && params.filter_sv_coverage) {
        mosdepth_cov_analysis_subworkflow(
            ch_input_bam_bai_bed,
            [[:], []]
        )
    }

    // Methylation calling with modkit (if enabled)
    ch_empty_bed = Channel.value([[:], []])

    if (params.methyl) {
        log.info "Running methylation calling with modkit"
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
    ================================================================================
                            STRUCTURAL VARIANT CALLING
    ================================================================================
    */
    
    if (params.sv) {
        // Run three different SV callers in parallel for comprehensive SV detection
        log.info "Running structural variant calling with Sniffles, CuteSV, and SVIM"
        // Sniffles: Population-scale SV calling
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

        /*
        ================================================================================
                                SV FILTERING BY COVERAGE
        ================================================================================
        */
        
        // Apply coverage-based filtering to SV calls if requested
        if (params.filter_sv_coverage) {
            log.info "Applying coverage-based filtering to SV calls"
            // Filter Sniffles results based on coverage depth and callable regions
            filterbycov_sniffles_subworkflow (
                sniffles_sv_subworkflow.out.vcf,
                mosdepth_subworkflow.out.summary_txt,
                mosdepth_cov_analysis_subworkflow.out.quantized_bed,
                params.chromosome_codes,
                params.min_read_support ?: 'auto',
                params.min_read_support_limit ?: 3,
                params.filter_pass ?: false
            )
            
            // Filter CuteSV results using same coverage criteria
            filterbycov_cutesv_subworkflow (
                cutesv_sv_subworkflow.out.vcf,
                mosdepth_subworkflow.out.summary_txt,
                mosdepth_cov_analysis_subworkflow.out.quantized_bed,
                params.chromosome_codes,
                params.min_read_support ?: 'auto',
                params.min_read_support_limit ?: 3,
                params.filter_pass ?: false
            )
            
            // Filter SVIM results using coverage information
            filterbycov_svim_subworkflow (
                svim_sv_subworkflow.out.vcf,
                mosdepth_subworkflow.out.summary_txt,
                mosdepth_cov_analysis_subworkflow.out.quantized_bed,
                params.chromosome_codes,
                params.min_read_support ?: 'auto',
                params.min_read_support_limit ?: 3,
                params.filter_pass ?: false
            )
            
            // Use filtered VCFs for downstream analysis
            ch_filtered_sniffles_vcf = filterbycov_sniffles_subworkflow.out.filtered_vcf
            ch_filtered_cutesv_vcf   = filterbycov_cutesv_subworkflow.out.filtered_vcf
            ch_filtered_svim_vcf     = filterbycov_svim_subworkflow.out.filtered_vcf
            
        } else {
            // Use original unfiltered VCFs
            ch_filtered_sniffles_vcf = sniffles_sv_subworkflow.out.vcf
            ch_filtered_cutesv_vcf   = cutesv_sv_subworkflow.out.vcf
            ch_filtered_svim_vcf     = svim_sv_subworkflow.out.vcf
        }

        /*
        ================================================================================
                            SV MERGING AND INTERSECTION FILTERING
        ================================================================================
        */
        
        // Merge SV calls from different callers using SURVIVOR or use individual results
        if (params.multi_caller_sv_filtering) {
            log.info "Combining SV calls from multiple callers and filtering using SURVIVOR"
            // Create chromosome sizes channel for bedtools intersect
            ch_chrom_sizes = params.chrom_sizes ? 
                Channel.value([[:], file(params.chrom_sizes)]) : 
                Channel.value([[:], []])

            // Prepare VCFs for SURVIVOR merging using the sample_name parameter
            ch_vcfs_for_merging = ch_filtered_sniffles_vcf
                .map { meta, vcf -> 
                    [params.sample_name, vcf, 'sniffles'] 
                }
                .mix(
                    ch_filtered_cutesv_vcf.map { meta, vcf -> 
                        [params.sample_name, vcf, 'cutesv'] 
                    }
                )
                .mix(
                    ch_filtered_svim_vcf.map { meta, vcf -> 
                        [params.sample_name, vcf, 'svim'] 
                    }
                )
                .groupTuple(by: 0)
                .map { sample_id, vcfs, callers ->
                    def meta = [id: sample_id]
                    [meta, vcfs]
                }

            // Run the multi-caller SV filtering subworkflow
            multi_caller_sv_filter_subworkflow(
                ch_vcfs_for_merging,           // VCFs for SURVIVOR merging
                ch_filtered_cutesv_vcf,        // Individual CuteSV VCFs
                ch_filtered_sniffles_vcf,      // Individual Sniffles VCFs  
                ch_filtered_svim_vcf,          // Individual SVIM VCFs
                ch_chrom_sizes,                // Chromosome sizes for bedtools
                params.max_distance_breakpoints,
                params.min_supporting_callers,
                params.account_for_type,
                params.account_for_sv_strands,
                params.estimate_distanced_by_sv_size,
                params.min_sv_size
            )

            // Use the subworkflow outputs for downstream analysis
            ch_concatenated_vcf = multi_caller_sv_filter_subworkflow.out.concatenated_vcf

        } else {  
            ch_concatenated_vcf = Channel.empty()
        }
        
    } else {
        // No SV calling - create empty channels
        ch_filtered_sniffles_vcf = Channel.empty()
        ch_filtered_cutesv_vcf   = Channel.empty()
        ch_filtered_svim_vcf     = Channel.empty()
        ch_concatenated_vcf      = Channel.empty()
    }

    /*
    ================================================================================
                            SINGLE NUCLEOTIDE VARIANT CALLING
    ================================================================================
    */
    
    if (params.snv) {
    log.info "Running SNV calling with Clair3"
    
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

    // Access outputs correctly
    ch_snv_vcf = clair3_snv_subworkflow.out.vcf
    ch_snv_tbi = clair3_snv_subworkflow.out.tbi

    if (params.filter_pass) {
        // Prepare input channels for BCFTOOLS_VIEW
        ch_vcf_index = ch_snv_vcf.join(ch_snv_tbi, by: 0)
        
        // Create empty channels for optional inputs
        ch_regions = Channel.value([])
        ch_targets = Channel.value([])  
        ch_samples = Channel.value([])
        
        // Call BCFTOOLS_VIEW with correct input structure
        bcftools_filter_clair3_subworkflow(
            ch_vcf_index,    // tuple val(meta), path(vcf), path(index)
            ch_regions,      // path(regions) - empty
            ch_targets,      // path(targets) - empty  
            ch_samples       // path(samples) - empty
        )
        
        // Use filtered results
        ch_final_snv_vcf = bcftools_filter_clair3_subworkflow.out.vcf
    } 
    else {
        // Use unfiltered results
        ch_final_snv_vcf = ch_snv_vcf
    }
    } 
    else {
    // No SNV calling - create empty channel
    ch_final_snv_vcf = Channel.empty()
    }
    

    // Run SNV calling with Clair3 and optionally DeepVariant
    if (params.run_deepvariant) {
    log.info "Running DeepVariant SNV calling"
    
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
        ch_vcf_index_dv = ch_dv_vcf.join(ch_dv_tbi, by: 0)
        
        // Create empty channels for optional inputs
        ch_regions = Channel.value([])
        ch_targets = Channel.value([])  
        ch_samples = Channel.value([])
        
        // Call BCFTOOLS_VIEW with correct input structure
        bcftools_filter_deepvariant_subworkflow(
            ch_vcf_index_dv, // tuple val(meta), path(vcf), path(index)
            ch_regions,      // path(regions) - empty
            ch_targets,      // path(targets) - empty  
            ch_samples       // path(samples) - empty
        )
        
        ch_final_dv_vcf = bcftools_filter_deepvariant_subworkflow.out.vcf
    } else {
        ch_final_dv_vcf = ch_dv_vcf
    }

    if (params.unify_snv_geneyx) {
        log.info "Unifying SNV VCFs from Clair3 and DeepVariant"
        if (!params.snv){
            error "Unifying SNV VCFs requires SNV calling to be enabled (--snv true)"
        }
        
        // Combine and concatenate VCF files from both callers
        combined_vcfs = ch_final_snv_vcf
            .join(ch_snv_tbi, by: 0)
            .join(ch_final_dv_vcf.join(ch_dv_tbi, by: 0), by: 0)
            .map { meta, clair3_vcf, clair3_tbi, dv_vcf, dv_tbi ->
                [
                    meta,
                    [clair3_vcf, dv_vcf],    // List of VCF files
                    [clair3_tbi, dv_tbi]     // List of corresponding TBI files
                ]
            }

        bcftools_concat_snv_subworkflow(combined_vcfs)
        ch_final_unified_snv = bcftools_concat_snv_subworkflow.out.vcf
    } else {
        ch_final_unified_snv = Channel.empty()
    }
    } else {
    ch_final_dv_vcf = Channel.empty()
    ch_final_unified_snv = Channel.empty()
    }
    
    /*
    ================================================================================
                                PHASING ANALYSIS
    ================================================================================
    */
    
    // Run phasing with LongPhase if enabled
    /*
================================================================================
                            PHASING ANALYSIS
================================================================================
*/

// Run phasing with LongPhase if enabled
if (params.phase && params.snv) {
    log.info "Running phasing analysis with LongPhase"
    
    if (params.sv && params.phase_with_sv) {
        // SV calling is enabled and want to use SVs for phasing
        
        // Determine which SV channel to use for phasing
        if (params.phase_with_multicalledSV) {
            ch_sv_for_phasing = ch_concatenated_vcf
        } else {
            log.info "Using individual SV caller VCFs for phasing"
            if (params.phase_sv_caller == 'sniffles') {
                log.info "Using Sniffles VCF for phasing"
                ch_sv_for_phasing = ch_filtered_sniffles_vcf
            } else if (params.phase_sv_caller == 'cutesv') {
                log.info "Using CuteSV VCF for phasing"
                ch_sv_for_phasing = ch_filtered_cutesv_vcf
            } else if (params.phase_sv_caller == 'svim') {
                log.info "Using SVIM VCF for phasing"
                ch_sv_for_phasing = ch_filtered_svim_vcf
            } else {
                error "Unknown SV caller: ${params.phase_sv_caller}. Supported: sniffles, cutesv, svim"
            }
        }

        // Phasing with both SNVs and SVs
        ch_longphase_input = ch_input_bam
            .join(ch_final_snv_vcf, by: 0)
            .join(ch_sv_for_phasing, by: 0)
            .map { meta, bam, bai, snv_vcf, sv_vcf -> 
                tuple(meta, bam, bai, snv_vcf, sv_vcf, []) 
            }

    } else {
        // Either SV calling is disabled OR we don't want to use SVs for phasing
        // Phasing with SNVs only
        ch_longphase_input = ch_input_bam
            .join(ch_final_snv_vcf, by: 0)
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
================================================================================
                        COPY NUMBER VARIANT CALLING
================================================================================
*/

if (params.cnv) {
    log.info "Running copy number variant calling"
    if (params.run_qdnaseq) {
        log.info "Running QDNAseq CNV calling"
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
        
    } else {
        // Spectre CNV calling - requires SNV data unless using test data
        if (params.use_test_data) {
            log.info "Running Spectre CNV calling in test mode with hardcoded parameters"
            // Test mode with hardcoded parameters
            ch_spectre_reference = Channel
                .fromPath(params.spectre_snv_vcf, checkIfExists: true)
                .map { vcf_file -> 
                    return [id: params.sample_name]  // Use sample_name parameter
                }
                .combine(Channel.fromPath(params.spectre_fasta_file, checkIfExists: true))
                .map { meta, fasta -> tuple(meta, fasta) }
            
            spectre_cnv_subworkflow(
                params.spectre_mosdepth,
                ch_spectre_reference,
                params.spectre_snv_vcf,
                params.spectre_metadata,
                params.spectre_blacklist
            )
            ch_spectre_vcf = spectre_cnv_subworkflow.out.vcf
        } else {
            log.info "Running Spectre CNV calling in production mode with Clair3 results"
            // Production mode using actual Clair3 pipeline results
            // Error handling: Check if SNV calling was enabled and produced results
            if (!params.snv) {
                error "Spectre CNV calling requires SNV calling to be enabled (--snv true) or use QDNAseq (--run_qdnaseq true)"
            }
            
            // Prepare reference channel for Spectre CNV calling
            ch_spectre_reference = ch_final_snv_vcf
                .map { meta, vcf_file -> 
                    return [id: params.sample_name]  // Use sample_name parameter
                }
                .combine(ch_fasta.map { meta, fasta -> fasta })
                .map { meta, fasta -> tuple(meta, fasta) }
            
            spectre_cnv_subworkflow(
                mosdepth_subworkflow.out.regions_bed,  
                ch_spectre_reference,
                ch_final_snv_vcf,      // Use Clair3 results 
                params.spectre_metadata,
                params.spectre_blacklist
            )
            ch_spectre_vcf = spectre_cnv_subworkflow.out.vcf
        }

        // Round decimal places in Spectre output
        round_dp_spectre_subworkflow(ch_spectre_vcf)
        ch_cnv_vcf = round_dp_spectre_subworkflow.out.vcf
    }
    } 
    else {
    ch_cnv_vcf = Channel.empty()
    }

/*
================================================================================
                        SHORT TANDEM REPEAT ANALYSIS
================================================================================
*/

if (params.str) {
    log.info "Running short tandem repeat analysis with STRaGLeR"
    straglr_str_subworkflow(
        ch_input_bam,
        ch_fasta,
        params.str_bed_file
    )
    ch_str_vcf = straglr_str_subworkflow.out.vcf
} else {
    ch_str_vcf = Channel.empty()
}


/*
================================================================================
                        VCF UNIFICATION (ALL THREE REQUIRED)
================================================================================
*/

if (params.unify_vcfs && params.sv && params.cnv && params.str) {
    log.info "Running VCF unification for SV, CNV, and STR calls"
    
    // Prepare individual VCF channels with consistent metadata
    ch_sv_for_unify = params.multi_caller_sv_filtering ? 
        ch_concatenated_vcf.map { meta, vcf -> [[id: params.sample_name], vcf] } :
        ch_filtered_sniffles_vcf.map { meta, vcf -> [[id: params.sample_name], vcf] }
    
    ch_cnv_for_unify = ch_cnv_vcf.map { meta, vcf -> [[id: params.sample_name], vcf] }
    
    ch_str_for_unify = ch_str_vcf.map { meta, vcf -> [[id: params.sample_name], vcf] }
    
    // Run VCF unification subworkflow
    unify_vcf_subworkflow(
        ch_sv_for_unify,                    // SV VCFs
        ch_cnv_for_unify,                   // CNV VCFs
        ch_str_for_unify,                   // STR VCFs
        params.modify_str_calls             // Modify repeat calls flag
    )
    
    ch_unified_vcf = unify_vcf_subworkflow.out.unified_vcf
    
} else if (params.unify_vcfs) {
    // Error if unify_vcfs is requested but not all three variant types are enabled
    error "VCF unification requires all three variant calling types to be enabled: --sv true --cnv true --str true"
} else {
    ch_unified_vcf = Channel.empty()
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