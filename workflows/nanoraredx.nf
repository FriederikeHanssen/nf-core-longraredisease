/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { bam_merge_workflow } from '../subworkflows/local/bam_merge.nf'
include { bam_to_fastq_workflow } from '../subworkflows/local/bam_to_fastq.nf'
include {ALIGN_MINIMAP2} from '../subworkflows/local/align_minimap2.nf'
include { sniffles_workflow } from '../subworkflows/local/sniffles_svc.nf'
include { cutesv_workflow    } from '../subworkflows/local/cutesv_svc.nf'
include { svim_workflow      } from '../subworkflows/local/svim_svc.nf'
include {mosdepth_workflow} from '../subworkflows/local/mosdepth.nf'
include { survivor_merge_workflow } from '../subworkflows/local/survivor_merge.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow nanoraredx {

    // ch_input_bam = Channel
    //     .fromPath(params.bam_file, checkIfExists: true)
    //     .map { bam ->
    //         def bai = file("${bam}.bai")
    //         def meta = [ id: bam.getBaseName().replaceFirst(/\.bam$/, '') ]
    //         tuple(meta, bam, bai)
    //     }

    // Input parameters - reference genome and tandem repeat file
    ch_fasta = Channel
        .fromPath(params.fasta_file, checkIfExists: true)
        .map { fasta -> tuple([id: "ref"], fasta) }

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
    bam_merge_workflow(ch_bam_files, [[:], []], [[:], []])
    
    // Convert merged unmapped BAMs to FASTQ using samtools fastq with methylation tags
    bam_to_fastq_workflow(bam_merge_workflow.out.bam)
    
    // Optional: View the FASTQ outputs
    // bam_to_fastq_workflow.out.other.view { meta, fastq_files ->
    //     "Sample ${meta.id}: Generated FASTQ files from merged BAM: ${fastq_files}"
    // }

   // Align fastq, sort and index  
   ALIGN_MINIMAP2(
        ch_fasta,
        bam_to_fastq_workflow.out.other
    )

    // Collect the sorted BAM and BAI files from the alignment step
    ch_input_bam = ALIGN_MINIMAP2.out.ch_sorted_bam
    .join(ALIGN_MINIMAP2.out.ch_sorted_bai, by: 0)
    .map { meta, bam, bai -> tuple(meta, bam, bai) }


    // Parallel execution of SV calling subworkflows using original input BAMs
    sniffles_workflow(
        ch_input_bam,
        ch_fasta,
        ch_trf,
        params.vcf_output,
        params.snf_output
    )

    cutesv_workflow(
        ch_input_bam,
        ch_fasta
    )

    svim_workflow(
        ch_input_bam,
        ch_fasta
    )

    // Collect all SV VCF files from the subworkflows

    //  ch_vcfs_grouped = sniffles_workflow.out.vcf
    //  .mix(cutesv_workflow.out.vcf)
    //  .mix(svim_workflow.out.vcf)
    //  .groupTuple()
    //  .map { meta, vcfs -> tuple(meta, vcfs) }

    // Merge SV VCF files using Survivor

    // survivor_merge_workflow(
    //     ch_vcfs_grouped,
    //     params.max_distance_breakpoints,
    //     params.min_supporting_callers,
    //     params.account_for_type,
    //     params.account_for_sv_strands,
    //     params.estimate_distanced_by_sv_size,
    //     params.min_sv_size
    // )

    // Collect the merged VCF file
    // survivor_merge_workflow.out.vcf.view { meta, vcf -> 
    //     "Merged VCF for ${meta.id}: ${vcf}"
    // }

    // Calculate coverage using mosdepth
    ch_input_mosdepth = ALIGN_MINIMAP2.out.ch_sorted_bam
    .join(ALIGN_MINIMAP2.out.ch_sorted_bai, by: 0)
    .map { meta, bam, bai ->
        def bed = params.bed_file ? file(params.bed_file) : []
        tuple(meta, bam, bai, bed)
    }


     mosdepth_workflow(
        ch_input_mosdepth,
        [[:], []]
     )



}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/