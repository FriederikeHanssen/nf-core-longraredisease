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
include { samtools_merge_bam_subworkflow     } from '../subworkflows/local/samtools_merge_bam.nf'
include { samtools_bam_to_fastq_subworkflow  } from '../subworkflows/local/samtools_bam_to_fastq.nf'
include { minimap2_align_subworkflow         } from '../subworkflows/local/minimap2_align.nf'
include { nanoplot_subworkflow           } from '../subworkflows/local/nanoplot.nf'

// Methylatio calling 

include { modkit_methyl_subworkflow } from '../subworkflows/local/modkit_methyl.nf'
include { modkit_cpg_subworkflow    } from '../subworkflows/local/modkit_cpg.nf'

// Coverage analysis subworkflows
include { mosdepth_subworkflow         } from '../subworkflows/local/mosdepth.nf'
include { mosdepth_cov_analysis_subworkflow  } from '../subworkflows/local/mosdepth_cov_analysis.nf'

// Structural variant calling subworkflows
include { sniffles_sv_subworkflow            } from '../subworkflows/local/sniffles_sv.nf'
include { cutesv_sv_subworkflow              } from '../subworkflows/local/cutesv_sv.nf'
include { svim_sv_subworkflow                } from '../subworkflows/local/svim_sv.nf'

// SV filtering subworkflows - coverage-based filtering
include { filterbycov_svim_subworkflow     } from '../subworkflows/local/filterbycov_svim.nf'
include { filterbycov_sniffles_subworkflow } from '../subworkflows/local/filterbycov_sniffles.nf'
include { filterbycov_cutesv_subworkflow   } from '../subworkflows/local/filterbycov_cutesv.nf'

// SV merging and intersection filtering subworkflows
include { multi_caller_sv_filter_subworkflow } from '../subworkflows/local/multi_caller_sv_filter.nf'

// SNV calling and processing subworkflows
include { clair3_snv_subworkflow             } from '../subworkflows/local/clair3_snv.nf'
include { deepvariant_snv_subworkflow        } from '../subworkflows/local/deepvariant_snv.nf'
include { bcftools_concat_snv_subworkflow    } from '../subworkflows/local/bcftools_concat_snv.nf'

// Phasing subworkflow
include { longphase_subworkflow              } from '../subworkflows/local/longphase.nf'

// CNV calling subworkflows
include { spectre_cnv_subworkflow            } from '../subworkflows/local/spectre_cnv.nf'
include { round_dp_spectre_subworkflow       } from '../subworkflows/local/round_dp_spectre.nf'
include { qdnaseq_cnv_subworkflow            } from '../subworkflows/local/qdnaseq_cnv.nf'

// STR analysis subworkflow
include { straglr_str_subworkflow            } from '../subworkflows/local/straglr_str.nf'

// VCF processing subworkflows
include { unify_vcf_subworkflow              } from '../subworkflows/local/unify_vcf.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow nanoraredx {

    /*
    ================================================================================
                                REFERENCE FILES SETUP
    ================================================================================
    */
    
    // Reference genome FASTA file
    ch_fasta = Channel
        .fromPath(params.fasta_file, checkIfExists: true)
        .map { fasta -> tuple([id: "ref"], fasta) }

    // Reference genome index file
    ch_fai = Channel
        .fromPath("${params.fasta_file}.fai", checkIfExists: true)
        .map { fai -> tuple([id: "ref_fai"], fai) }

    ch_fasta_fai = Channel
    .fromPath(params.fasta_file, checkIfExists: true)
    .map { fasta -> 
        def fai = file("${fasta}.fai")
        tuple([id: "ref"], fasta, fai)
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
    
    // Collect and group BAM files by sample ID for merging
    ch_bam_files = Channel
        .fromPath("${params.bam_dir}/*.bam")
        .map { bam ->
            def sample_id = bam.name.split('_')[0]  // Extract sample ID from filename
            return [ [id: sample_id], bam ]
        }
        .groupTuple()

    // Merge multiple BAM files per sample into single BAM
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

    /*
    ================================================================================
                                COVERAGE ANALYSIS
    ================================================================================
    */
    
    // Prepare input channel with BAM, BAI, and optional BED file for coverage analysis
    ch_input_bam_bai_bed = minimap2_align_subworkflow.out.ch_sorted_bam
        .join(minimap2_align_subworkflow.out.ch_sorted_bai, by: 0)
        .map { meta, bam, bai ->
            def bed = params.bed_file ? file(params.bed_file) : []
            tuple(meta, bam, bai, bed)
        }

    // Calculate depth statistics using mosdepth (required for filtering if enabled)
    if (params.sv && params.filter_sv_coverage) {
        mosdepth_subworkflow(
            ch_input_bam_bai_bed,
            [[:], []]
        )

        // Perform coverage analysis and quantization for filtering
        mosdepth_cov_analysis_subworkflow(
            ch_input_bam_bai_bed,
            [[:], []]
        )
    }

    // Prepare input for nanoplot 
    ch_fastq_nanoplot = samtools_bam_to_fastq_subworkflow.out.other
    .map { meta, fastq_file ->
        tuple(meta, fastq_file) }

    
    // Prepare simplified BAM input channel for variant calling and methylation calling
    
    ch_input_bam = minimap2_align_subworkflow.out.ch_sorted_bam
        .join(minimap2_align_subworkflow.out.ch_sorted_bai, by: 0)
        .map { meta, bam, bai -> tuple(meta, bam, bai) }


    // Run nanoplot 

    if (params.qc) {
    nanoplot_subworkflow(
        ch_fastq_nanoplot
    )
    }


    // Methylation calling with modkit (if enabled)

    ch_empty_bed = Channel.value([[:], []])

    if (params.methyl){
    modkit_methyl_subworkflow(
        ch_input_bam,
        ch_fasta_fai,
        ch_empty_bed  // Empty channel for optional BED file
    )

    modkit_cpg_subworkflow(
        ch_input_bam,
        ch_fasta_fai,
        ch_empty_bed // Empty channel for optional BED file
    ) }

    /*
    ================================================================================
                            STRUCTURAL VARIANT CALLING
    ================================================================================
    */
    
    if (params.sv) {
        // Run three different SV callers in parallel for comprehensive SV detection
        
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
            
            // Filter Sniffles results based on coverage depth and callable regions
            filterbycov_sniffles_subworkflow (
                sniffles_sv_subworkflow.out.vcf,
                mosdepth_subworkflow.out.summary_txt,
                mosdepth_cov_analysis_subworkflow.out.quantized_bed,
                params.chromosome_codes,
                params.min_read_support ?: 'auto',
                params.min_read_support_limit ?: 2
            )
            
            // Filter CuteSV results using same coverage criteria
            filterbycov_cutesv_subworkflow (
                cutesv_sv_subworkflow.out.vcf,
                mosdepth_subworkflow.out.summary_txt,
                mosdepth_cov_analysis_subworkflow.out.quantized_bed,
                params.chromosome_codes,
                params.min_read_support ?: 'auto',
                params.min_read_support_limit ?: 2
            )
            
            // Filter SVIM results using coverage information
            filterbycov_svim_subworkflow (
                svim_sv_subworkflow.out.vcf,
                mosdepth_subworkflow.out.summary_txt,
                mosdepth_cov_analysis_subworkflow.out.quantized_bed,
                params.chromosome_codes,
                params.min_read_support ?: 'auto',
                params.min_read_support_limit ?: 3
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
            
            // Create chromosome sizes channel for bedtools intersect
            ch_chrom_sizes = params.chrom_sizes ? 
                Channel.value([[:], file(params.chrom_sizes)]) : 
                Channel.value([[:], []])

            // Prepare VCFs for SURVIVOR merging by extracting base sample IDs
            ch_vcfs_for_merging = ch_filtered_sniffles_vcf
                .map { meta, vcf -> 
                    def sample_id = meta.id.split('_')[0]
                    [sample_id, vcf, 'sniffles'] 
                }
                .mix(
                    ch_filtered_cutesv_vcf.map { meta, vcf -> 
                        def sample_id = meta.id.split('_')[0]
                        [sample_id, vcf, 'cutesv'] 
                    }
                )
                .mix(
                    ch_filtered_svim_vcf.map { meta, vcf -> 
                        def sample_id = meta.id.split('_')[0]
                        [sample_id, vcf, 'svim'] 
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
            // Use individual caller results without merging/intersection filtering
            ch_final_sniffles_vcf = ch_filtered_sniffles_vcf
            ch_final_cutesv_vcf   = ch_filtered_cutesv_vcf
            ch_final_svim_vcf     = ch_filtered_svim_vcf
            ch_concatenated_vcf   = Channel.empty()
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

        // Run SNV calling with Clair3 and optionally DeepVariant
        if (params.use_deepvariant) {
            
            // Run both Clair3 and DeepVariant for comprehensive SNV detection
            deepvariant_snv_subworkflow(
                ch_input_bam_bai_bed,
                ch_fasta,
                ch_fai,
                [[:], []],
                [[:], []]
            )

            clair3_snv_subworkflow(
                ch_input_bam_clair3,
                ch_fasta,
                ch_fai
            )

            // Combine and concatenate VCF files from both callers
            combined_vcfs = clair3_snv_subworkflow.out.vcf
                .join(clair3_snv_subworkflow.out.tbi, by: 0)
                .join(
                    deepvariant_snv_subworkflow.out.vcf
                        .join(deepvariant_snv_subworkflow.out.vcf_tbi, by: 0),
                    by: 0
                )
                .map { meta, clair3_vcf, clair3_tbi, dv_vcf, dv_tbi ->
                    [
                        meta,
                        [clair3_vcf, dv_vcf],    // List of VCF files
                        [clair3_tbi, dv_tbi]     // List of corresponding TBI files
                    ]
                }

            results_snv = bcftools_concat_snv_subworkflow(combined_vcfs)
            
        } else {
            // Run only Clair3 SNV calling
            results_snv = clair3_snv_subworkflow(
                ch_input_bam_clair3,
                ch_fasta,
                ch_fai
            )
        }
    } else {
        // No SNV calling - create empty channel
        results_snv = [vcf: Channel.empty()]
    }

    /*
    ================================================================================
                                PHASING ANALYSIS
    ================================================================================
    */
    
    // Run phasing with LongPhase if enabled
    if (params.phase && params.snv) {

        // Determine which SV channel to use for phasing
        if (params.phase_with_sv && params.sv) {

            if (params.multi_caller_sv_filtering) {
                // Use merged and filtered SV VCF from multi-caller consensus
                ch_sv_for_phasing = ch_concatenated_vcf
            } else {
                // Use parameter to dynamically select the channel
                if (params.phase_sv_caller_unify == 'sniffles') {
                    ch_sv_for_phasing = ch_filtered_sniffles_vcf
                } else if (params.phase_sv_caller_unify == 'cutesv') {
                    ch_sv_for_phasing = ch_filtered_cutesv_vcf
                } else if (params.phase_sv_caller_unify == 'svim') {
                    ch_sv_for_phasing = ch_filtered_svim_vcf
                } else {
                    error "Unknown SV caller: ${params.phase_sv_caller_unify}. Supported: sniffles, cutesv, svim"
                }
            }

            // Phasing with both SNVs and SVs
            ch_longphase_input = ch_input_bam
                .join(results_snv.vcf, by: 0)
                .join(ch_sv_for_phasing, by: 0)
                .map { meta, bam, bai, snv_vcf, sv_vcf -> 
                    tuple(meta, bam, bai, snv_vcf, sv_vcf, []) 
                }

        } else {
            // Phasing with SNVs only
            ch_longphase_input = ch_input_bam
                .join(results_snv.vcf, by: 0)
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
        
        if (params.use_qdnaseq) {
            // CNV calling with QDNAseq (suitable for targeted/exome data)
            results_cnv = qdnaseq_cnv_subworkflow(
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
            // CNV calling with Spectre (requires whole genome data)
            ch_spectre_reference = Channel
                .fromPath(params.spectre_snv_vcf, checkIfExists: true)
                .map { vcf_file -> 
                    def sample_id = vcf_file.baseName.split('_')[0]
                    return [id: sample_id]
                }
                .combine(Channel.fromPath(params.spectre_fasta_file, checkIfExists: true))
                .map { meta, fasta -> tuple(meta, fasta) }
            
            results_spectre = spectre_cnv_subworkflow(
                params.spectre_mosdepth,
                ch_spectre_reference,
                params.spectre_snv_vcf,
                params.spectre_metadata,
                params.spectre_blacklist
            )

            // Optional: Round decimal places in Spectre output
            results_cnv = round_dp_spectre_subworkflow(results_spectre.vcf)
        }
        
    } else {
        results_cnv = [vcf: Channel.empty()]
    }

    /*
    ================================================================================
                            SHORT TANDEM REPEAT ANALYSIS
    ================================================================================
    */
    
    if (params.str) {
        results_str = straglr_str_subworkflow(
            ch_input_bam,
            ch_fasta,
            params.str_bed_file
        )
    } else {
        results_str = [vcf: Channel.empty()]
    }


/*
================================================================================
                            VCF UNIFICATION
================================================================================
*/



if (params.unify_vcfs && (params.sv || params.cnv || params.str)) {
    
    // Prepare SV VCF channel (only if SV analysis is enabled)
    ch_sv_for_unify = params.sv ? 
        (params.multi_caller_sv_filtering ? 
            ch_concatenated_vcf.map { meta, vcf -> 
                def sample_id = meta.id.split('_')[0]
                def new_meta = [id: sample_id]
                [new_meta, vcf]
            } :
            ch_filtered_sniffles_vcf.map { meta, vcf ->
                def sample_id = meta.id.split('_')[0]
                def new_meta = [id: sample_id]
                [new_meta, vcf]
            }
        ) :
        Channel.empty()
    
    // Prepare CNV VCF channel (only if CNV analysis is enabled)
    ch_cnv_for_unify = params.cnv ? 
        results_cnv.vcf.map { meta, vcf ->
            def sample_id = meta.id.split('_')[0]
            def new_meta = [id: sample_id]
            [new_meta, vcf]
        } :
        Channel.empty()
    
    // Prepare STR VCF channel (only if STR analysis is enabled)
    ch_str_for_unify = params.str ? 
        results_str.vcf.map { meta, vcf ->
            def sample_id = meta.id.split('_')[0]
            def new_meta = [id: sample_id]
            [new_meta, vcf]
        } :
        Channel.empty()
    
    // Get all unique sample IDs from enabled analyses
    ch_all_samples = Channel.empty()
    if (params.sv) {
        ch_all_samples = ch_all_samples.mix(
            ch_sv_for_unify.map { meta, vcf -> meta }
        )
    }
    if (params.cnv) {
        ch_all_samples = ch_all_samples.mix(
            ch_cnv_for_unify.map { meta, vcf -> meta }
        )
    }
    if (params.str) {
        ch_all_samples = ch_all_samples.mix(
            ch_str_for_unify.map { meta, vcf -> meta }
        )
    }
    
    // Get unique sample metas
    ch_sample_metas = ch_all_samples.unique()
    
    // Create channels with empty files for missing VCF types
    ch_empty_file = Channel.value([[:], []])
    
    // Join all channels by sample ID, using remainder: true and providing empty defaults
    ch_unified_input = ch_sample_metas
        .join(ch_sv_for_unify, by: 0, remainder: true)
        .join(ch_cnv_for_unify, by: 0, remainder: true)
        .join(ch_str_for_unify, by: 0, remainder: true)
        .map { meta, sv_vcf, cnv_vcf, str_vcf ->
            [
                [meta, sv_vcf ?: []],      // Empty list if no SV VCF
                [meta, cnv_vcf ?: []],     // Empty list if no CNV VCF  
                [meta, str_vcf ?: []]      // Empty list if no STR VCF
            ]
        }
    
    // Split into separate channels
    ch_sv_final = ch_unified_input.map { sv, cnv, str -> sv }
    ch_cnv_final = ch_unified_input.map { sv, cnv, str -> cnv }
    ch_str_final = ch_unified_input.map { sv, cnv, str -> str }
    
    // Run VCF unification subworkflow
    unify_vcf_subworkflow(
        ch_sv_final,                        // SV VCFs (or empty list)
        ch_cnv_final,                       // CNV VCFs (or empty list)
        ch_str_final,                       // STR VCFs (or empty list)
        params.modify_str_calls             // Modify repeat calls flag
    )
    
    // Set final unified output
    ch_unified_vcf = unify_vcf_subworkflow.out.unified_vcf
    
} else {
    // No unification requested or no analysis types enabled
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