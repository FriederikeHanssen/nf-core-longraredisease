#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NANORAREDX: Comprehensive Nanopore Rare Disease Analysis Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Data preprocessing subworkflows
include { BAM_STATS_SAMTOOLS                 } from '../subworkflows/nf-core/bam_stats_samtools/main.nf'
include { SAMTOOLS_INDEX                     } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FAIDX                     } from '../modules/nf-core/samtools/faidx/main.nf'
include { bam2fastq_subworkflow              } from '../subworkflows/local/bam2fastq.nf'
include { minimap2_align_bam_subworkflow         } from '../subworkflows/local/minimap2_align_bam.nf'
include { minimap2_align_fastq_subworkflow         } from '../subworkflows/local/minimap2_align_fastq.nf'
include { CAT_FASTQ                          } from '../modules/nf-core/cat/fastq/main.nf'
include { NANOPLOT as NANOPLOT_QC            } from '../modules/nf-core/nanoplot/main'

// Methylation calling 
include { methyl_subworkflow                 } from '../subworkflows/local/methylation.nf'

// Coverage analysis subworkflows
include { mosdepth_subworkflow               } from '../subworkflows/local/mosdepth.nf'

// Structural variant calling subworkflows
include { sv_subworkflow                      } from '../subworkflows/local/sv.nf'
include { SVANNA_PRIORITIZE                   } from '../modules/local/SvAnna/main.nf'

// SV merging and intersection filtering subworkflows

include { GUNZIP as GUNZIP_SNIFFLES           } from '../modules/nf-core/gunzip/main.nf'
include { GUNZIP as GUNZIP_CUTESV           } from '../modules/nf-core/gunzip/main.nf'
include { GUNZIP as GUNZIP_SVIM           } from '../modules/nf-core/gunzip/main.nf'
include { consensuSV_subworkflow              } from '../subworkflows/local/consensuSV.nf'


// SNV calling and processing subworkflows
include { snv_subworkflow                     } from '../subworkflows/local/snv.nf'
include { merge_snv_subworkflow               } from '../subworkflows/local/merge_snv.nf'

// Phasing subworkflow
include { longphase_subworkflow                } from '../subworkflows/local/longphase.nf'

// CNV calling subworkflows
include { cnv_subworkflow                      } from '../subworkflows/local/cnv.nf'
include { cnv_qdnaseq_subworkflow               } from '../subworkflows/local/cnv_qdnaseq.nf'

// STR analysis subworkflow
include { str_subworkflow              } from '../subworkflows/local/str.nf'

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

if (params.filter_sv && !params.sv) {
    error "SV coverage filtering requires SV calling to be enabled (--sv true)"
}

if (params.cnv && !params.run_qdnaseq && !params.use_test_data && !params.snv) {
    error "Spectre CNV calling requires SNV calling to be enabled (--snv true) unless using test data"
}

if (params.consensuSV && !params.sv) {
    error "Multi-caller SV filtering requires SV calling to be enabled (--sv true)"
}

// ADD: Check for conflicting alignment options
if (params.align_with_bam && params.align_with_fastq) {
    error "Cannot use both --align_with_bam and --align_with_fastq simultaneously. Choose one alignment method."
}

// Input validation based on alignment parameter
if (params.align_with_bam) {
    // When align_with_bam=true, expect unaligned BAMs in bam_dir
    if (!params.bam_dir) {
        error "When --align_with_bam is true, --bam_dir must be provided with unaligned BAM files"
    }
    if (!params.sample_name) {
        error "When --align_with_bam is true, --sample_name must be provided"
    }
} else if (params.align_with_fastq) {
    // When align_with_fastq=true, expect FASTQ files in fastq_dir
    if (!params.fastq_dir) {
        error "When --align_with_fastq is true, --fastq_dir must be provided with FASTQ files"
    }
    if (!params.sample_name) {
        error "When --align_with_fastq is true, --sample_name must be provided"
    }
} else {
    // When neither alignment option is true, expect aligned BAM file
    if (!params.aligned_bam) {
        error "When alignment is disabled, --aligned_bam must be provided with the path to aligned BAM file"
    }
    // Extract sample name from aligned BAM filename if not provided
    if (!params.sample_name) {
        def aligned_bam_file = file(params.aligned_bam)
        def extracted_sample_name = aligned_bam_file.name
            .replaceAll(/\.sorted\.bam$/, '')
            .replaceAll(/\.bam$/, '')
        
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
    SAMTOOLS_FAIDX(ch_fasta, true)
    ch_fai = SAMTOOLS_FAIDX.out.fai

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
    if (params.align_with_fastq) {
        /*
        ================================================================================
                            FASTQ ALIGNMENT WORKFLOW
        ================================================================================
        */
        
        log.info "Sample name: ${params.sample_name}"
        
        // Collect FASTQ files
        ch_fastq_files = Channel
            .fromPath("${params.fastq_dir}/*.fastq.gz", checkIfExists: true)
            .collect()  // Collect all files into a single list
            .map { fastq_list ->
                // Create single sample with all FASTQ files
                def meta = [id: params.sample_name]
                return [meta, fastq_list]
            }

        // Align FASTQ reads to reference genome using minimap2
        minimap2_align_fastq_subworkflow(
            ch_fasta,
            ch_fastq_files
        )

        // Set final aligned BAM channels from minimap2 output
        ch_final_sorted_bam = minimap2_align_fastq_subworkflow.out.bam
        ch_final_sorted_bai = minimap2_align_fastq_subworkflow.out.index

        // Prepare input for nanoplot from FASTQ
        CAT_FASTQ(
            ch_fastq_files.map { meta, fastq_list ->
        [meta + [single_end: true], fastq_list]  // Add single_end: true inline
        }
        )
        ch_nanoplot = CAT_FASTQ.out.reads
        
    }

    else if (params.align_with_bam) {
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

        // Convert BAM to FASTQ
        bam2fastq_subworkflow(
            ch_bam_files, 
            [[:], []], 
            [[:], []]
        )

        // Align FASTQ reads to reference genome using minimap2
        minimap2_align_bam_subworkflow(
            ch_fasta,
            bam2fastq_subworkflow.out.other
        )

        // Set final aligned BAM channels from minimap2 output
        ch_final_sorted_bam = minimap2_align_bam_subworkflow.out.ch_sorted_bam
        ch_final_sorted_bai = minimap2_align_bam_subworkflow.out.ch_sorted_bai

        // Prepare input for nanoplot from FASTQ
        ch_nanoplot = bam2fastq_subworkflow.out.other
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
        // Use pre-aligned BAM file (assume BAI exists)
        ch_final_sorted_bam = Channel
        .fromPath(params.aligned_bam, checkIfExists: true)
        .map { bam ->
        def meta = [id: params.sample_name]
        return [meta, bam]
        }
        
        ch_final_sorted_bai = ch_final_sorted_bam.map { meta, bam ->
        return [meta, file("${bam}.bai")]
        }
        // For nanoplot, we'll skip it (no FASTQ available)
        ch_nanoplot = ch_final_sorted_bam
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
    if (params.qc) {
        // Both workflows produce ch_fastq_nanoplot
        NANOPLOT_QC(
            ch_nanoplot
        )
    }
    
    // Run mosdepth when needed
    if (params.sv || params.cnv || params.generate_coverage) {
        mosdepth_subworkflow(
            ch_input_bam_bai_bed,
            [[:], []]
        )
    }

    // Methylation calling with modkit (if enabled)
    ch_empty_bed = Channel.value([[:], []])



    if (params.methyl) {
        if (params.align_with_bam) {
        // Use workflow-generated BAM for methylation analysis
            ch_methyl_input = ch_input_bam
            
        
        } else {
            // For align_with_fastq or no alignment: use methylated BAM from path
            if (!params.methyl_bam) {
                error "When --methyl is enabled without --align_with_bam, --methyl_bam must be provided as fastq files do not have methylation info"
            }
        
            // Create methylation BAM channel following the same logic as ch_input_bam
            ch_methyl_input = Channel
            .fromPath(params.methyl_bam, checkIfExists: true)
            .map { bam ->
            def meta = [id: params.sample_name]
            return [meta, bam, file("${bam}.bai")] 
            }
            }
            
            methyl_subworkflow(
            ch_methyl_input,
            ch_fasta_fai,
            ch_empty_bed  // Empty channel for optional BED file
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
        
        // Run SV calling subworkflow
        sv_subworkflow(
        ch_input_bam,
        ch_fasta,
        ch_trf,
        params.vcf_output,
        params.snf_output,
        params.primary_sv_caller,
        params.filter_sv,
        mosdepth_subworkflow.out.summary_txt,
        mosdepth_subworkflow.out.quantized_bed,
        params.chromosome_codes ?: 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY',
        params.min_read_support ?: 'auto',
        params.min_read_support_limit ?: 3,
        params.filter_pass ?: false
        )
        

        // Extract VCF from the sv_gz_tbi channel for unify_vcf_subworkflow
            ch_sv_vcf = sv_subworkflow.out.primary_vcf_gz
                

        /*
        ================================================================================
                            MULTI-CALLER FILTERING AND CONSENSUS
        ================================================================================
        */
        
        if (params.consensuSV) {
            // Prepare VCFs for SURVIVOR merging - direct from subworkflow outputs
            ch_vcfs_for_merging = GUNZIP_SNIFFLES(sv_subworkflow.out.sniffles_vcf_gz).gunzip
                .map { meta, vcf -> 
                    [params.sample_name ?: meta.id, vcf, 'sniffles'] 
                }
                .mix(
                    GUNZIP_CUTESV(sv_subworkflow.out.cutesv_vcf_gz).gunzip
                        .filter { meta, vcf -> vcf && vcf.exists() }
                        .map { meta, vcf -> 
                            [params.sample_name ?: meta.id, vcf, 'cutesv'] 
                        }
                )
                .mix(
                    GUNZIP_SVIM(sv_subworkflow.out.svim_vcf_gz).gunzip
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


                // Group by meta to ensure we're merging files from the same sample
                ch_merge_input = sv_subworkflow.out.sniffles_vcf_gz
                .join(sv_subworkflow.out.sniffles_tbi, by: 0)
                .join(sv_subworkflow.out.svim_vcf_gz, by: 0)
                .join(sv_subworkflow.out.svim_tbi, by: 0)
                .join(sv_subworkflow.out.cutesv_vcf_gz, by: 0)
                .join(sv_subworkflow.out.cutesv_tbi, by: 0)
                .map { meta, sniffles_vcf, sniffles_tbi, svim_vcf, svim_tbi, cutesv_vcf, cutesv_tbi ->
                tuple(meta, 
                  [sniffles_vcf, svim_vcf, cutesv_vcf], 
                  [sniffles_tbi, svim_tbi, cutesv_tbi])
                  }

                // Run multi-caller filtering
                consensuSV_subworkflow(
                ch_vcfs_for_merging,
                ch_merge_input
                )

            // Extract VCF from the gz_tbi channel for unify_vcf_subworkflow
            // In your consensuSV section, normalize the meta before emitting
            ch_sv_vcf = consensuSV_subworkflow.out.vcf
            .map { meta, vcf_gz -> 
             def clean_meta = [id: meta.id]
             tuple(clean_meta, vcf_gz) 
             }
            
        }
    } 
     
     else {
        // Create empty channel when SV calling is disabled
        ch_sv_vcf = Channel.empty()
    }

    SVANNA_PRIORITIZE(
    ch_sv_vcf,           // tuple val(meta), path(vcf)
    params.svanna_db,    // path(data_directory) 
    params.hpo_terms     // val(hpo_terms)
    )

/*
================================================================================
                        SINGLE NUCLEOTIDE VARIANT CALLING
================================================================================
*/

    if (params.snv) {
        // Prepare input for SNV calling
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
        
        // Run SNV calling
        snv_subworkflow(
            ch_input_bam_clair3,
            ch_fasta,
            ch_fai,
            params.deepvariant,
            ch_input_bam_bai_bed,
            params.filter_pass
        )

        ch_snv_vcf = snv_subworkflow.out.clair3_vcf
        ch_snv_tbi = snv_subworkflow.out.clair3_tbi

        if (params.merge_snv && params.deepvariant) {
            combined_vcfs = ch_snv_vcf
                .join(ch_snv_tbi, by: 0)
                .join(
                    snv_subworkflow.out.deepvariant_vcf
                        .join(snv_subworkflow.out.deepvariant_tbi, by: 0), 
                    by: 0
                )
                .map { meta, clair3_vcf, clair3_tbi, deepvariant_vcf, deepvariant_tbi ->
                    [
                        meta,
                        [clair3_vcf, deepvariant_vcf],    // List of VCF files
                        [clair3_tbi, deepvariant_tbi]     // List of corresponding TBI files
                    ]
                }
        
            // Merge SNV VCFs
            merge_snv_subworkflow(combined_vcfs)
        }
    } else {
        // Create empty channels when SNV calling is disabled
        ch_snv_vcf = Channel.empty()
        ch_snv_tbi = Channel.empty()
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
                .join(ch_snv_vcf, by: 0)
                .join(ch_sv_vcf, by: 0)
                .map { meta, bam, bai, snv_vcf, sv_vcf -> 
                    tuple(meta, bam, bai, snv_vcf, sv_vcf, []) 
                }
        } else {
            // Phasing with SNVs only
            ch_longphase_input = ch_input_bam
                .join(ch_snv_vcf, by: 0)
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
            // Test mode with hardcoded parameters
            ch_spectre_test_reference = Channel
                .fromPath(params.spectre_test_clair3_vcf, checkIfExists: true)
                .map { vcf_file -> 
                    return [id: params.sample_name]  // Use sample_name parameter
                }
                .combine(Channel.fromPath(params.spectre_test_fasta_file, checkIfExists: true))
                .map { meta, fasta -> tuple(meta, fasta) }
            
            cnv_subworkflow(
                params.spectre_test_mosdepth,
                ch_spectre_test_reference,
                params.spectre_test_clair3_vcf,
                params.spectre_metadata,
                params.spectre_blacklist
            )
            ch_cnv_vcf = cnv_subworkflow.out.vcf
        } else {
            // Extract just the BED file from mosdepth output
            ch_spectre_mosdepth_bed = mosdepth_subworkflow.out.regions_bed.map { meta, bed -> bed }
            
            // Prepare reference channel for Spectre CNV calling
            ch_spectre_reference = ch_snv_vcf
                .map { meta, vcf_file -> 
                    [id: params.sample_name]  // Create new meta with sample_name
                }
                .combine(Channel.fromPath(params.fasta_file, checkIfExists: true))
                .map { meta, fasta -> tuple(meta, fasta) }
            
            // Extract just the VCF file from the ch_snv_vcf channel
            ch_spectre_clair3_vcf = ch_snv_vcf.map { meta, vcf -> vcf }
    
            cnv_subworkflow(
                ch_spectre_mosdepth_bed,
                ch_spectre_reference, 
                ch_spectre_clair3_vcf,              
                params.spectre_metadata,
                params.spectre_blacklist
            )

            ch_cnv_vcf = cnv_subworkflow.out.vcf
        }
    } else {
        // Create empty channel when CNV calling is disabled
        ch_cnv_vcf = Channel.empty()
    }
    
    if (params.cnv && params.run_qdnaseq) {
        // QDNAseq doesn't need mosdepth or SNV data
        cnv_qdnaseq_subworkflow(
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
        str_subworkflow(
            ch_input_bam,
            ch_fasta,
            params.str_bed_file
        )
        ch_str_vcf = str_subworkflow.out.vcf
    } else {
        ch_str_vcf = Channel.empty()
    }
  
/*
================================================================================
                             VCF UNIFICATION 
================================================================================
*/

   
        if (params.unify_geneyx) {
            
            unify_vcf_subworkflow(
            params.sv ? ch_sv_vcf : Channel.value([[:], []]),
            params.cnv ? ch_cnv_vcf : Channel.value([[:], []]),
            params.str ? ch_str_vcf : Channel.value([[:], []]),
            params.modify_str_calls ?: false
        )
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


