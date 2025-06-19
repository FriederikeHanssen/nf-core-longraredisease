/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { samtools_merge_bam_subworkflow } from '../subworkflows/local/samtools_merge_bam.nf'
include { samtools_bam_to_fastq_subworkflow  } from '../subworkflows/local/samtools_bam_to_fastq.nf'
include { minimap2_align_subworkflow } from '../subworkflows/local/minimap2_align.nf'
include { mosdepth_cnv_depth_subworkflow } from '../subworkflows/local/mosdepth_cnv_depth.nf'
//include { samtools_bam_downsample_subworkflow } from '../subworkflows/local/samtools_bam_downsample.nf' // optional
// include { nanoplot_subworkflow } from '../subworkflows/local/nanoplot.nf'

include { sniffles_sv_subworkflow } from '../subworkflows/local/sniffles_sv.nf'
include { cutesv_sv_subworkflow    } from '../subworkflows/local/cutesv_sv.nf'
include { svim_sv_subworkflow      } from '../subworkflows/local/svim_sv.nf'
include { survivor_merge_sv_subworkflow } from '../subworkflows/local/survivor_merge_sv.nf'
// include { filterbycov_sv_subworkflow } from '../subworkflows/local/filtercalls_sv.nf'
// include { filter_vcf_insdeldup_subworkflow } from '../subworkflows/local/filter_vcf_insdeldup.nf'
// include {savanna_sv_annotation_subworkflow} from '../subworkflows/local/savanna_sv_annotation.nf'

// include { clair3_model_selection_subworkflow } from '../subworkflows/local/clair3_model_selection.nf'
include { clair3_snv_subworkflow } from '../subworkflows/local/clair3_snv.nf'
include { deepvariant_snv_subworkflow } from '../subworkflows/local/deepvariant_snv.nf'
include { bcftools_concat_snv_subworkflow } from '../subworkflows/local/bcftools_concat_snv.nf'
// include {longphase_snp_subworkflow} from '../subworkflows/local/longphase_snp.nf'

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
    

// Optional: View the FASTQ outputs
    // bam_to_fastq_workflow.out.other.view { meta, fastq_files ->
    //     "Sample ${meta.id}: Generated FASTQ files from merged BAM: ${fastq_files}"
    // }


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

// Downsample bam files if specified 


// Collect the sorted BAM and BAI files from the alignment step depending if downsampling is enabled or not
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

// Boolean - filter calls by coverage
// Boolean - filter calls by type
// Merge SV calls from different callers using Survivor
// Collect all SV VCF files from the subworkflows

    ch_vcfs_grouped = sniffles_sv_subworkflow.out.vcf
    .mix(svim_sv_subworkflow.out.vcf)
    .mix(cutesv_sv_subworkflow.out.vcf)
    .groupTuple()
    .map { meta, vcfs -> tuple(meta, vcfs) }

    // Merge SV VCF files using Survivor

    survivor_merge_sv_subworkflow(
        ch_vcfs_grouped,
        params.max_distance_breakpoints,
        params.min_supporting_callers,
        params.account_for_type,
        params.account_for_sv_strands,
        params.estimate_distanced_by_sv_size,
        params.min_sv_size
    )



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

// Choose your own Clair3 model selection otherwise run teh standard before running Clair3 SNV calling
    

// Run deepvariant SNV calling on the aligned BAM files - choose deepvariant/ clair3 or both
// Run SNV calling using Clair3 or DeepVariant depending on --snv value

// if (params.snv) {
//     if (params.use_deepvariant) {
//         deepvariant_snv_subworkflow(
//             ch_input_bam_bai_bed,
//             ch_fasta,
//             ch_fai,
//             [[:], []],
//             [[:], []])

//         clair3_snv_subworkflow(
//             ch_input_bam_clair3,
//             ch_fasta,
//             ch_fai)

//         // Combine VCF files with their TBI indices
//         // Assuming clair3_snv_subworkflow also has .vcf and .tbi outputs
//         combined_vcfs = clair3_snv_subworkflow.out.vcf
//             .join(clair3_snv_subworkflow.out.tbi, by: 0)  // Join Clair3 VCF with its TBI
//             .join(
//                 deepvariant_snv_subworkflow.out.vcf
//                     .join(deepvariant_snv_subworkflow.out.vcf_tbi, by: 0),  // Join DeepVariant VCF with its TBI
//                 by: 0
//             )
//             .map { meta, clair3_vcf, clair3_tbi, dv_vcf, dv_tbi ->
//                 [
//                     meta,
//                     [clair3_vcf, dv_vcf],  // List of VCF files
//                     [clair3_tbi, dv_tbi]   // List of corresponding TBI files
//                 ]
//             }

//         results_snv = bcftools_concat_snv_subworkflow(combined_vcfs)
        
//         } 
        
//         else {
//         // Only run Clair3 SNV calling
//         results_snv = clair3_snv_subworkflow(
//             ch_input_bam_clair3,
//             ch_fasta,
//             ch_fai)
//             }
//             } 

//             else {
//                 results_snv = Channel.empty()
//                 }

    clair3_snv_subworkflow(
            ch_input_bam_clair3,
            ch_fasta,
            ch_fai)

   deepvariant_snv_subworkflow(
        ch_input_bam_bai_bed,
        ch_fasta,
        ch_fai,
        [[:], []],
        [[:], []])

// Run phasing with LongPhase on the aligned BAM files - optional if enabled  
   
// Run CNV calling with Spectre or QDNASeq

// Input parameters - reference genome and tandem repeat file

    ch_spectre_reference = Channel
    .fromPath(params.spectre_snv_vcf, checkIfExists: true)
    .map { vcf_file -> 
        def sample_id = vcf_file.baseName.split('_')[0]
        return [id: sample_id]
    }
    .combine(Channel.fromPath(params.spectre_fasta_file, checkIfExists: true))
    .map { meta, fasta -> tuple(meta, fasta) }

// If spectre is selected then run RoundDP to remove NAs for compatibility with Geneyx

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
        // cnv calling with spectre

        } else {
            results_cnv =  spectre_cnv_subworkflow(
        params.spectre_mosdepth,
        params.spectre_bin_size,
        ch_spectre_reference,
        params.spectre_snv_vcf,
        params.spectre_metadata,
        params.spectre_blacklist
    )
        }
        }

    else {
        results_cnv = Channel.empty()
        }



        

//  Run STR calling with STRaglr

    ch_input_bam_str = Channel.fromPath(params.str_bam, checkIfExists: true)
    .map { bam ->
        def bai = file("${bam}.bai")  // Assumes BAI is named same as BAM + .bai
        def meta = [id: bam.getBaseName().replaceFirst(/\.bam$/, '')]
        tuple(meta, bam, bai)
    }

    if (params.str) {
        // use haplotagged bam from snp() as input to str()
        bam_channel_str = ch_input_bam_str

        results_str = straglr_str_subworkflow(
          bam_channel_str,
          ch_fasta,
          params.str_bed_file
        )
        
        } else {

        results_str = Channel.empty()
        }

// Methylation calling with Modkit

// Unify VCFs from different workflows

    //     if (!params.snp && !params.sv && !params.mod && !params.cnv && !params.str) {
    //     log.error (colors.red + "No work to be done! Choose one or more workflows to run from [--snp, --sv, --cnv, --str, --mod]" + colors.reset)
    //     can_start = false
    // }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/