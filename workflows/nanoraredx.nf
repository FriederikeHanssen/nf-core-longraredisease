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

// Import nf-schema function
include { samplesheetToList } from 'plugin/nf-schema'

// Data preprocessing subworkflows
include { BAM_STATS_SAMTOOLS                 } from '../subworkflows/nf-core/bam_stats_samtools/main.nf'
include { SAMTOOLS_INDEX                     } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FAIDX                     } from '../modules/nf-core/samtools/faidx/main.nf'
include { bam2fastq_subworkflow              } from '../subworkflows/local/bam2fastq.nf'
include { minimap2_align_subworkflow     } from '../subworkflows/local/minimap2_align.nf'
include { CAT_FASTQ                          } from '../modules/nf-core/cat/fastq/main.nf'
include { NANOPLOT as NANOPLOT_QC            } from '../modules/nf-core/nanoplot/main'

// Methylation calling 
include { methyl_subworkflow                 } from '../subworkflows/local/methylation.nf'

// Coverage analysis subworkflows
include { mosdepth_subworkflow               } from '../subworkflows/local/mosdepth.nf'

// Structural variant calling subworkflows
include { sv_subworkflow                     } from '../subworkflows/local/sv.nf'
include { SVANNA_PRIORITIZE                  } from '../modules/local/SvAnna/main.nf'

// SV merging and intersection filtering subworkflows
include { GUNZIP as GUNZIP_SNIFFLES          } from '../modules/nf-core/gunzip/main.nf'
include { GUNZIP as GUNZIP_CUTESV            } from '../modules/nf-core/gunzip/main.nf'
include { GUNZIP as GUNZIP_SVIM              } from '../modules/nf-core/gunzip/main.nf'
include { consensuSV_subworkflow             } from '../subworkflows/local/consensuSV.nf'

// SNV calling and processing subworkflows
include { snv_subworkflow                    } from '../subworkflows/local/snv.nf'
include { merge_snv_subworkflow              } from '../subworkflows/local/merge_snv.nf'

// Phasing subworkflow
include { longphase_subworkflow              } from '../subworkflows/local/longphase.nf'

// CNV calling subworkflows
include { cnv_subworkflow                    } from '../subworkflows/local/cnv.nf'
include { cnv_qdnaseq_subworkflow            } from '../subworkflows/local/cnv_qdnaseq.nf'

// STR analysis subworkflow
include { str_subworkflow                    } from '../subworkflows/local/str.nf'

// VCF processing subworkflows
include { unify_vcf_subworkflow              } from '../subworkflows/local/unify_vcf.nf'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_nanoraredx_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow nanoraredx {

    main:
    
    // Convert samplesheet to list and create channel using nf-schema
    def samplesheet_data = samplesheetToList(params.input, "assets/schema_input.json")
    
    ch_samplesheet = Channel.fromList(samplesheet_data)
        .map { row ->
            // Handle the ArrayList structure from nf-schema
            if (row instanceof List) {
                def meta_map = row[0]
                def sample_id = meta_map.id ?: meta_map.toString()
                def meta = [id: sample_id]
                def data = [
                    bam_dir: row[1] ?: null,
                    fastq_dir: row[2] ?: null, 
                    aligned_bam: row[3] ?: null,
                    methyl_bam: row[4] ?: null,
                    hpo_terms: row[5] ?: null
                ]
                return [meta, data]
            } else {
                error "Unexpected row type: ${row.getClass()}"
            }
        }

/*
=======================================================================================
                                REFERENCE FILES SETUP
=======================================================================================
*/
    // Initialize versions channel
    ch_versions = Channel.empty()

    ch_fasta = Channel
        .fromPath(params.fasta_file, checkIfExists: true)
        .map { fasta -> tuple([id: "ref"], fasta) }
        .first()

    // Generate FAI index
    SAMTOOLS_FAIDX(ch_fasta, true)
    ch_fai = SAMTOOLS_FAIDX.out.fai
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    // Create combined FASTA+FAI channel by joining
    ch_fasta_fai = ch_fasta
        .join(ch_fai, by: 0)
        .map { meta, fasta, fai -> tuple(meta, fasta, fai) }
        .first()
    
    // Tandem repeat file for Sniffles (only if SV calling is enabled)
    if (params.sv) {
        ch_trf = Channel
            .fromPath(params.sniffles_tandem_file, checkIfExists: true)
            .map { bed -> tuple([id: "trf"], bed) }
            .first()
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
        // Collect FASTQ files
        ch_samplesheet.view()
        ch_fastq_files = ch_samplesheet
            .map { meta, data ->
                def fastq_dir = file(data.fastq_dir)
                def fastq_files = fastq_dir.listFiles().findAll { 
                    it.name.endsWith('.fastq.gz') || it.name.endsWith('.fq.gz') 
                }
                return [meta, fastq_files]
            }

        // Align FASTQ reads to reference genome using minimap2
        minimap2_align_subworkflow(
            ch_fasta,
            ch_fastq_files
        )
        ch_versions = ch_versions.mix(minimap2_align_subworkflow.out.versions)


        // Set final aligned BAM channels from minimap2 output
        ch_final_sorted_bam = minimap2_align_subworkflow.out.bam
        ch_final_sorted_bai = minimap2_align_subworkflow.out.bai

        // Prepare input for nanoplot from FASTQ
        CAT_FASTQ(
            ch_fastq_files.map { meta, fastq_list ->
                [meta + [single_end: true], fastq_list]
            }
        )
        ch_nanoplot = CAT_FASTQ.out.reads
        ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)
    }

    else if (params.align_with_bam) {
        /*
        ================================================================================
                            ALIGNMENT WORKFLOW (UNALIGNED INPUT)
        ================================================================================
        */
        
        // Collect unaligned BAM files
        ch_bam_files = ch_samplesheet
            .map { meta, data ->
                def bam_pattern = "${data.bam_dir}/*.bam"
                def bam_files = file(bam_pattern)
                
                // Ensure bam_files is always a list
                def bam_list = bam_files instanceof List ? bam_files : [bam_files]
                
                if (bam_list.isEmpty()) {
                    error "No BAM files found for sample ${meta.id} in directory: ${data.bam_dir}"
                }
                
                return [meta, bam_list]
            }

        // Convert BAM to FASTQ
        bam2fastq_subworkflow(
            ch_bam_files, 
            [[:], []], 
            [[:], []]
        )

        ch_versions = ch_versions.mix(bam2fastq_subworkflow.out.versions)
        // Align FASTQ reads to reference genome using minimap2
        minimap2_align_subworkflow(
            ch_fasta,
            bam2fastq_subworkflow.out.other
        )
        ch_versions = ch_versions.mix(minimap2_align_subworkflow.out.versions)
        // Set final aligned BAM channels from minimap2 output
        ch_final_sorted_bam = minimap2_align_subworkflow.out.bam
        ch_final_sorted_bai = minimap2_align_subworkflow.out.bai

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
        
        // For aligned BAM input
        ch_aligned_input = ch_samplesheet
            .map { meta, data ->
                def bam_file = file(data.aligned_bam, checkIfExists: true)
                def bai_file = file("${data.aligned_bam}.bai", checkIfExists: true)
                return [meta, bam_file, bai_file]
            }

        // Use this single channel for all downstream processes
        ch_final_sorted_bam = ch_aligned_input.map { meta, bam, bai -> [meta, bam] }
        ch_final_sorted_bai = ch_aligned_input.map { meta, bam, bai -> [meta, bai] }
        
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
        ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)
    }
    
    // Run nanoplot (only if we have FASTQ data from alignment workflow)
    if (params.qc) {
        NANOPLOT_QC(
            ch_nanoplot
        )
        ch_versions = ch_versions.mix(NANOPLOT_QC.out.versions)
    }
    
    // Run mosdepth when needed
    if (params.sv || params.cnv || params.generate_coverage) {
        mosdepth_subworkflow(
            ch_input_bam_bai_bed,
            [[:], []]
        )
        ch_versions = ch_versions.mix(mosdepth_subworkflow.out.versions)
    }


    if (params.methyl) {
        if (params.align_with_bam) {
            // Use workflow-generated BAM for methylation analysis
            ch_methyl_input = ch_input_bam
        } else {
            // For align_with_fastq or no alignment: use methylated BAM from path
            ch_methyl_input = ch_samplesheet
                .map { meta, data ->
                    if (!data.methyl_bam) {
                        error "When --methyl is enabled without --align_with_bam, methyl_bam must be provided in samplesheet for sample ${meta.id}"
                    }
                    def bam_file = file(data.methyl_bam, checkIfExists: true)
                    def bai_file = file("${data.methyl_bam}.bai", checkIfExists: true)
                    return [meta, bam_file, bai_file]
                }
        }
        
        methyl_subworkflow(
            ch_methyl_input,
            ch_fasta_fai,
            [[:], []]
        )
        ch_versions = ch_versions.mix(methyl_subworkflow.out.versions)
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
            params.filter_pass_sv ?: false
        )

        ch_versions = ch_versions.mix(sv_subworkflow.out.versions)

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
                    [meta.id, vcf, 'sniffles'] 
                }
                .mix(
                    GUNZIP_CUTESV(sv_subworkflow.out.cutesv_vcf_gz).gunzip
                        .filter { meta, vcf -> vcf && vcf.exists() }
                        .map { meta, vcf -> 
                            [meta.id, vcf, 'cutesv'] 
                        }
                )
                .mix(
                    GUNZIP_SVIM(sv_subworkflow.out.svim_vcf_gz).gunzip
                        .filter { meta, vcf -> vcf && vcf.exists() }
                        .map { meta, vcf -> 
                            [meta.id, vcf, 'svim'] 
                        }
                )
                .groupTuple(by: 0)
                .map { sample_id, vcfs, callers ->
                    def meta = [id: sample_id, callers: callers]
                    [meta, vcfs]
                }

            ch_versions = ch_versions.mix(GUNZIP_SNIFFLES.out.versions)
            ch_versions = ch_versions.mix(GUNZIP_CUTESV.out.versions)
            ch_versions = ch_versions.mix(GUNZIP_SVIM.out.versions)

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
                ch_merge_input,
                params.use_survivor_bed
            )
            ch_versions = ch_versions.mix(consensuSV_subworkflow.out.versions)


            ch_sv_vcf = consensuSV_subworkflow.out.vcf
                .map { meta, vcf_gz -> 
                    def clean_meta = [id: meta.id]
                    tuple(clean_meta, vcf_gz) 
                }
        }

    if (params.annotate_sv) {
    
         ch_hpo_terms = ch_samplesheet.map { meta, data -> [meta, data.hpo_terms] }

         ch_sv_vcf_with_hpo = ch_sv_vcf
            .join(ch_hpo_terms, by: 0)
    
         ch_svanna_db = Channel
         .fromPath(params.svanna_db, checkIfExists: true)
         .first()

         SVANNA_PRIORITIZE(
                 ch_sv_vcf_with_hpo.map { meta, vcf, hpo_terms -> [meta, vcf] },
                 ch_svanna_db,
                 ch_sv_vcf_with_hpo.map { meta, vcf, hpo_terms -> hpo_terms }
             )
         ch_versions = ch_versions.mix(SVANNA_PRIORITIZE.out.versions) 
    
    }

    } else {
        ch_sv_vcf = Channel.empty()
    }
/*
================================================================================
                        SINGLE NUCLEOTIDE VARIANT CALLING
================================================================================
*/

    if (params.snv) {
        // Prepare input for SNV calling
        ch_input_bam_clair3 = ch_input_bam.map { meta, bam, bai ->
            tuple(
                meta,
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
            params.filter_pass_snv
        )

        ch_versions = ch_versions.mix(snv_subworkflow.out.versions)

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
                        [clair3_vcf, deepvariant_vcf],
                        [clair3_tbi, deepvariant_tbi]
                    ]
                }
        
            // Merge SNV VCFs
            merge_snv_subworkflow(combined_vcfs)
            ch_versions = ch_versions.mix(merge_snv_subworkflow.out.versions)
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
         ch_versions = ch_versions.mix(longphase_subworkflow.out.versions)
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
            ch_spectre_test_reference = ch_samplesheet
            .map { meta, data -> meta.id }  // Extract sample ID from samplesheet
            .combine(Channel.fromPath(params.spectre_test_clair3_vcf, checkIfExists: true))
            .combine(Channel.fromPath(params.spectre_test_fasta_file, checkIfExists: true))
            .map { sample_id, vcf_file, fasta -> 
            def meta = [id: sample_id]
            tuple(meta, fasta)
            }
            
            cnv_subworkflow(
                params.spectre_test_mosdepth,
                ch_spectre_test_reference,
                params.spectre_test_clair3_vcf,
                params.spectre_metadata,
                params.spectre_blacklist
            )
            ch_cnv_vcf = cnv_subworkflow.out.vcf
            ch_versions = ch_versions.mix(cnv_subworkflow.out.versions)
        } 
        
        else {

        ch_combined = ch_snv_vcf
        .join(mosdepth_subworkflow.out.regions_bed, by: 0)
        // Result: [meta, vcf_file, bed_file]
  
        // Transform for cnv_subworkflow - assuming it expects separate channels
        ch_spectre_bed = ch_combined.map { meta, vcf, bed -> bed }
        ch_spectre_vcf = ch_combined.map { meta, vcf, bed -> vcf }

        ch_spectre_reference = ch_samplesheet
        .map { meta, data -> meta.id }  
        .join(
            ch_combined.map { meta, vcf, bed -> [meta.id, vcf] },  
            by: 0
        )  // Combine with VCF
        .combine(Channel.fromPath(params.fasta_file, checkIfExists: true))  
        .map { sample_id, vcf_file, fasta -> 
            def meta = [id: sample_id]
            tuple(meta, fasta)
        }
        
        cnv_subworkflow(
        ch_spectre_bed,
        ch_spectre_reference,
        ch_spectre_vcf,
        params.spectre_metadata,
        params.spectre_blacklist
        )

        ch_cnv_vcf = cnv_subworkflow.out.vcf
        ch_versions = ch_versions.mix(cnv_subworkflow.out.versions)
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
        ch_versions = ch_versions.mix(cnv_qdnaseq_subworkflow.out.versions)
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
        ch_versions = ch_versions.mix(str_subworkflow.out.versions)
    } else {
        ch_str_vcf = Channel.empty()
    }
  
/*
================================================================================
                             VCF UNIFICATION 
================================================================================
*/

    if (params.unify_geneyx) {

        ch_combined = ch_sv_vcf
        .join(ch_cnv_vcf, by: 0, remainder: true)
        .join(ch_str_vcf, by: 0, remainder: true)

        
        unify_vcf_subworkflow(
        ch_combined.map { meta, sv, cnv, str -> [meta, sv] },
        ch_combined.map { meta, sv, cnv, str -> [meta, cnv ?: []] },
        ch_combined.map { meta, sv, cnv, str -> [meta, str ?: []] },
        params.modify_str_calls ?: false

    )
    ch_versions = ch_versions.mix(unify_vcf_subworkflow.out.versions)
    

    //  unify_vcf_subworkflow(
        //     params.sv ? ch_sv_vcf : Channel.value([[:], []]),
        //     params.cnv ? ch_cnv_vcf : Channel.value([[:], []]),
        //     params.str ? ch_str_vcf : Channel.value([[:], []]),
        //     params.modify_str_calls ?: false
        // )

    }

    softwareVersionsToYAML(ch_versions)
    .collectFile(
        storeDir: "${params.outdir}/pipeline_info",
        name: 'nf_core_nanoraredx_software_versions.yml',
        sort: true,
        newLine: true
    ).set { ch_collated_versions }


emit:
    versions = ch_collated_versions


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/