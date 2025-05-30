/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { bam_merge_workflow } from '../subworkflows/local/bam_merge.nf'
include { bam_index_workflow } from '../subworkflows/local/bam_index.nf'
include { bam_to_fastq_workflow } from '../subworkflows/local/bam_to_fastq.nf'
include { sniffles_workflow } from '../subworkflows/local/sniffles_svc.nf'
include { cutesv_workflow    } from '../subworkflows/local/cutesv_svc.nf'
include { svim_workflow      } from '../subworkflows/local/svim_svc.nf'
include {ALIGN_MINIMAP2} from '../subworkflows/local/align_minimap2.nf'

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

    ch_fasta = Channel
        .fromPath(params.fasta_file, checkIfExists: true)
        .map { fasta -> tuple([id: "ref"], fasta) }
    
    // ch_fasta_mm = Channel
    //     .fromPath(params.fasta_mm_file, checkIfExists: true)
    //     .map { fasta_mm -> tuple([id: "ref_mm"], fasta_mm) }

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
    

    // Debug
    // ch_bam_files.view { meta, bams -> 
    //     "Sample: ${meta.id}, BAM files: ${bams}" 
    // }

    // Run BAM merge workflow
    bam_merge_workflow(ch_bam_files, [[:], []], [[:], []])
    
    // Convert merged unmapped BAMs to FASTQ using samtools fastq with methylation tags
    bam_to_fastq_workflow(bam_merge_workflow.out.bam)
    
   

    // Optional: View the FASTQ outputs
    bam_to_fastq_workflow.out.other.view { meta, fastq_files ->
        "Sample ${meta.id}: Generated FASTQ files from merged BAM: ${fastq_files}"
    }

   // Align fastq, sort and index  
   ALIGN_MINIMAP2(
        ch_fasta,
        bam_to_fastq_workflow.out.other
    )

    ch_input_bam = ALIGN_MINIMAP2.out.ch_sorted_bam
    .join(ALIGN_MINIMAP2.out.ch_sorted_bai, by: 0)
    .map { meta, bam, bai -> tuple(meta, bam, bai) }

        // ch_input_bam = Channel
    //     .fromPath(params.bam_file, checkIfExists: true)
    //     .map { bam ->
    //         def bai = file("${bam}.bai")
    //         def meta = [ id: bam.getBaseName().replaceFirst(/\.bam$/, '') ]
    //         tuple(meta, bam, bai)
    //     }

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
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/