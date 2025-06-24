/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { samtools_merge_bam_subworkflow } from '../subworkflows/local/samtools_merge_bam.nf'
include { samtools_bam_to_fastq_subworkflow  } from '../subworkflows/local/samtools_bam_to_fastq.nf'
include { minimap2_align_subworkflow } from '../subworkflows/local/minimap2_align.nf'
include { mosdepth_cnv_depth_subworkflow } from '../subworkflows/local/mosdepth_cnv_depth.nf'
// include { nanoplot_subworkflow } from '../subworkflows/local/nanoplot.nf'

include { sniffles_sv_subworkflow } from '../subworkflows/local/sniffles_sv.nf'
include { cutesv_sv_subworkflow    } from '../subworkflows/local/cutesv_sv.nf'
include { svim_sv_subworkflow      } from '../subworkflows/local/svim_sv.nf'
include { survivor_merge_sv_subworkflow } from '../subworkflows/local/survivor_merge_sv.nf'
include { filterbycov_sv_svim } from '../subworkflows/local/filtersv_svim.nf'
include { filterbycov_sv_sniffles } from '../subworkflows/local/filtersv_sniffles.nf'
include { filterbycov_sv_cutesv } from '../subworkflows/local/filtersv_cutesv.nf'
// include {savanna_sv_annotation_subworkflow} from '../subworkflows/local/savanna_sv_annotation.nf'

include { clair3_snv_subworkflow } from '../subworkflows/local/clair3_snv.nf'
include { deepvariant_snv_subworkflow } from '../subworkflows/local/deepvariant_snv.nf'
include { bcftools_concat_snv_subworkflow } from '../subworkflows/local/bcftools_concat_snv.nf'
include {longphase_subworkflow} from '../subworkflows/local/longphase.nf'

include { spectre_cnv_subworkflow } from '../subworkflows/local/spectre_cnv.nf'
// include {rounddp_spectre_str_subworkflow} from '../subworkflows/local/rounddp_spectre_str.nf'
include { qdnaseq_cnv_subworkflow } from '../subworkflows/local/qdnaseq_cnv'

include { straglr_str_subworkflow } from '../subworkflows/local/straglr_str.nf'

// include { modkit_mc_subworkflow } from '../subworkflows/local/modkit_mc.nf'
// include { meow_mc_subworkflow } from '../subworkflows/local/meow_mc.nf'
//  include { unifyvcf_sv_snv_str_subworkflow } from '../subworkflows/local/unifyvcf_sv_snv_str.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow nanoraredx {

//  Collect reference files
    ch_fasta = Channel
        .fromPath(params.fasta_file, checkIfExists: true)
        .map { fasta -> tuple([id: "ref"], fasta) }

    ch_fai = Channel
        .fromPath("${params.fasta_file}.fai", checkIfExists: true)
        .map { fai -> tuple([id: "ref_fai"], fai) }

    ch_trf = Channel
        .fromPath(params.tandem_file, checkIfExists: true)
        .map { bed -> tuple([id: "trf"], bed) }

// Collect all BAMs for merging
    ch_bam_files = Channel.fromPath("${params.bam_dir}/*.bam")
        .map { bam ->
            def sample_id = bam.name.split('_')[0]
            return [ [id: sample_id], bam ]
        }
        .groupTuple()

// Run BAM merge workflow
    samtools_merge_bam_subworkflow(ch_bam_files, [[:], []], [[:], []])

// Convert merged unmapped BAMs to FASTQ using samtools fastq with methylation tags
    samtools_bam_to_fastq_subworkflow(samtools_merge_bam_subworkflow.out.bam)

// Align fastq, sort and index  
    minimap2_align_subworkflow(
        ch_fasta,
        samtools_bam_to_fastq_subworkflow.out.other
    )

// Calculate coverage/cnv using mosdepth
    ch_input_bam_bai_bed = minimap2_align_subworkflow.out.ch_sorted_bam
        .join(minimap2_align_subworkflow.out.ch_sorted_bai, by: 0)
        .map { meta, bam, bai ->
            def bed = params.bed_file ? file(params.bed_file) : []
            tuple(meta, bam, bai, bed)
        }

    mosdepth_cnv_depth_subworkflow (
        ch_input_bam_bai_bed,
        [[:], []]
    )

    // Collect the sorted BAM and BAI files from the alignment step
    ch_input_bam = minimap2_align_subworkflow.out.ch_sorted_bam
        .join(minimap2_align_subworkflow.out.ch_sorted_bai, by: 0)
        .map { meta, bam, bai -> tuple(meta, bam, bai) }

// Parallel execution of SV calling subworkflows using original input BAMs
    sniffles_sv_subworkflow(
        ch_input_bam,
        ch_fasta,
        ch_trf,
        params.vcf_output,
        params.snf_output
    )

    cutesv_sv_subworkflow(
        ch_input_bam,
        ch_fasta
    )

    svim_sv_subworkflow(
        ch_input_bam,
        ch_fasta
    )

    if (params.filter_sv_calls) {
        // Filter individual SV caller results for all three callers
        
        // Filter Sniffles results
        filterbycov_sv_sniffles(
            sniffles_sv_subworkflow.out.vcf,
            mosdepth_cnv_depth_subworkflow.out.summary_txt,
            mosdepth_cnv_depth_subworkflow.out.quantized_bed,
            params.chromosome_codes?.split(',') ?: ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrMT'],
            params.min_read_support ?: 'auto',
            params.min_read_support_limit ?: 2
        )
        
        // Filter CuteSV results
        filterbycov_sv_cutesv(
            cutesv_sv_subworkflow.out.vcf,
            mosdepth_cnv_depth_subworkflow.out.summary_txt,
            mosdepth_cnv_depth_subworkflow.out.quantized_bed,
            params.chromosome_codes?.split(',') ?: ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrMT'],
            params.min_read_support ?: 'auto',
            params.min_read_support_limit ?: 2
        )
        
        // Filter SVIM results
        filterbycov_sv_svim(
            svim_sv_subworkflow.out.vcf,
            mosdepth_cnv_depth_subworkflow.out.summary_txt,
            mosdepth_cnv_depth_subworkflow.out.quantized_bed,
            params.chromosome_codes?.split(',') ?: ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrMT'],
            params.min_read_support ?: 'auto',
            params.min_read_support_limit ?: 2
        )
        
        // Use filtered VCFs for downstream analysis
        ch_filtered_sniffles_vcf = filterbycov_sv_sniffles.out.filtered_vcf
        ch_filtered_cutesv_vcf = filterbycov_sv_cutesv.out.filtered_vcf
        ch_filtered_svim_vcf = filterbycov_sv_svim.out.filtered_vcf
    } else {
        ch_filtered_sniffles_vcf = sniffles_sv_subworkflow.out.vcf
        ch_filtered_cutesv_vcf = cutesv_sv_subworkflow.out.vcf
        ch_filtered_svim_vcf = svim_sv_subworkflow.out.vcf
    }

// Merge SV calls from different callers using Survivor OR use individual caller results
    if (params.merge_sv_calls) {
    // Create a channel that groups by sample ID only
    ch_vcfs_for_merging = ch_filtered_sniffles_vcf
        .map { meta, vcf -> [meta.id, vcf, 'sniffles'] }
        .mix(
            ch_filtered_cutesv_vcf.map { meta, vcf -> [meta.id, vcf, 'cutesv'] }
        )
        .mix(
            ch_filtered_svim_vcf.map { meta, vcf -> [meta.id, vcf, 'svim'] }
        )
        .groupTuple(by: 0)
        .map { sample_id, vcfs, callers ->
            def meta = [id: sample_id]
            [meta, vcfs]
        }
    
    ch_vcfs_for_merging.view { "VCFs for merging: $it" }
    
    survivor_merge_sv_subworkflow(
        ch_vcfs_for_merging,
        params.max_distance_breakpoints,
        params.min_supporting_callers,
        params.account_for_type,
        params.account_for_sv_strands,
        params.estimate_distanced_by_sv_size,
        params.min_sv_size
    )
}

// Run Clair3 SNV calling on the aligned BAM files
    ch_input_bam_clair3 = ch_input_bam.map { meta, bam, bai ->
        def id = bam.getBaseName().replaceAll(/\.bam$/, '')  // Remove .bam extension
        def updated_meta = meta + [ id: id ]                 // Inject/override `id` into meta

        tuple(
            updated_meta,
            bam,
            bai,
            params.clair3_model,
            [],
            params.clair3_platform
        )
    }

// Run SNV calling using Clair3 or DeepVariant depending on --snv value
    if (params.use_deepvariant) {
        deepvariant_snv_subworkflow(
            ch_input_bam_bai_bed,
            ch_fasta,
            ch_fai,
            [[:], []],
            [[:], []])

        clair3_snv_subworkflow(
            ch_input_bam_clair3,
            ch_fasta,
            ch_fai)

        // Combine VCF files with their TBI indices
        combined_vcfs = clair3_snv_subworkflow.out.vcf
            .join(clair3_snv_subworkflow.out.tbi, by: 0)  // Join Clair3 VCF with its TBI
            .join(
                deepvariant_snv_subworkflow.out.vcf
                    .join(deepvariant_snv_subworkflow.out.vcf_tbi, by: 0),  // Join DeepVariant VCF with its TBI
                by: 0
            )
            .map { meta, clair3_vcf, clair3_tbi, dv_vcf, dv_tbi ->
                [
                    meta,
                    [clair3_vcf, dv_vcf],  // List of VCF files
                    [clair3_tbi, dv_tbi]   // List of corresponding TBI files
                ]
            }

        results_snv = bcftools_concat_snv_subworkflow(combined_vcfs)
    } else {
        // Only run Clair3 SNV calling
        results_snv = clair3_snv_subworkflow(
            ch_input_bam_clair3,
            ch_fasta,
            ch_fai)
    }

// Run phasing with LongPhase on the aligned BAM files - optional if enabled 
    if (params.phase & params.phase_with_sv) {
        // Determine which SV results to use based on merge_sv_calls parameter
        def sv_results_channel
        
        if (params.merge_sv_calls) {
            // Use merged SV results
            sv_results_channel = survivor_merge_sv_subworkflow.out.vcf
        } else {
            // Use individual SV caller results based on params.sv_caller (now using filtered versions)
            if (params.phase_sv_caller == 'sniffles') {
                sv_results_channel = ch_filtered_sniffles_vcf
            } else if (params.phase_sv_caller == 'cutesv') {
                sv_results_channel = ch_filtered_cutesv_vcf
            } else if (params.phase_sv_caller == 'svim') {
                sv_results_channel = ch_filtered_svim_vcf
            } else {
                error "Unknown SV caller: ${params.sv_caller}. Supported callers: sniffles, cutesv, svim"
            }
        }
        
        // Join BAM data with SNV results and SV results
        ch_longphase_input = ch_input_bam
            .join(results_snv.vcf, by: 0, remainder: true)
            .join(sv_results_channel, by: 0, remainder: true)
            .map { meta, bam, bai, snv_vcf, sv_vcf ->
                tuple(
                    meta,
                    bam,
                    bai,
                    snv_vcf,
                    sv_vcf,
                    []
                )
            }

        longphase_subworkflow(
            ch_longphase_input,
            ch_fasta,
            ch_fai
        )
    } else if (params.phase){
        ch_longphase_input = ch_input_bam
            .join(results_snv.vcf, by: 0, remainder: true)
            .map { meta, bam, bai, snv_vcf ->
                tuple(
                    meta,
                    bam,
                    bai,
                    snv_vcf,
                    [],
                    []
                )
            }

        longphase_subworkflow(
            ch_longphase_input,
            ch_fasta,
            ch_fai
        )

    }

// Run CNV calling with Spectre or QDNASeq
    ch_spectre_reference = Channel
        .fromPath(params.spectre_snv_vcf, checkIfExists: true)
        .map { vcf_file -> 
            def sample_id = vcf_file.baseName.split('_')[0]
            return [id: sample_id]
        }
        .combine(Channel.fromPath(params.spectre_fasta_file, checkIfExists: true))
        .map { meta, fasta -> tuple(meta, fasta) }

    if (params.cnv) {
        // cnv calling with qdnaseq
        if (params.use_qdnaseq) {
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
            // cnv calling with spectre - can not be run on test data or subset of data 
            // requires whole genome
            results_cnv = spectre_cnv_subworkflow(
                params.spectre_mosdepth,
                params.spectre_bin_size,
                ch_spectre_reference,
                params.spectre_snv_vcf,
                params.spectre_metadata,
                params.spectre_blacklist
            )
        }
    } else {
        results_cnv = Channel.empty()
    }

// Run STR calling with STRaglr
    if (params.str) {
        results_str = straglr_str_subworkflow(
            ch_input_bam,
            ch_fasta,
            params.str_bed_file
        )
    } else {
        results_str = Channel.empty()
    }

// Methylation calling with Modkit
// Unify VCFs from different workflows

    // if (!params.snp && !params.sv && !params.mod && !params.cnv && !params.str) {
    //     log.error (colors.red + "No work to be done! Choose one or more workflows to run from [--snp, --sv, --cnv, --str, --mod]" + colors.reset)
    //     can_start = false
    // }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/