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
include { survivor_filter_workflow } from '../subworkflows/local/survivor_filter_lowcov.nf'
include { clair3_snv_workflow } from '../subworkflows/local/clair3_snv.nf'
include { spectre_cnv_workflow } from '../subworkflows/local/spectre_cnv.nf'


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

// Run Clair3 SNV calling on the aligned BAM files
    ch_input_clair3 = ch_input_bam.map { meta, bam, bai ->
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

    
    clair3_snv_workflow(
        ch_input_clair3,
        ch_fasta,
        ch_fai
    )

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

// Debugging before calling SPECTRE
// clair3_snv_workflow.out.vcf.view { "VCF for SPECTRE: $it" }
// mosdepth_workflow.out.regions_bed.view { "Mosdepth files: $it" }
// ch_fasta.view { "Reference for SPECTRE: $it" }

// Collect all SV VCF files from the subworkflows

    ch_vcfs_grouped = sniffles_workflow.out.vcf
    .mix(svim_workflow.out.vcf)
    .mix(cutesv_workflow.out.vcf)
    .groupTuple()
    .map { meta, vcfs -> tuple(meta, vcfs) }

    // Merge SV VCF files using Survivor

    survivor_merge_workflow(
        ch_vcfs_grouped,
        params.max_distance_breakpoints,
        params.min_supporting_callers,
        params.account_for_type,
        params.account_for_sv_strands,
        params.estimate_distanced_by_sv_size,
        params.min_sv_size
    )


// Create conditional channels
// ch_for_filtering = params.run_survivor_filter ? 
//     survivor_merge_workflow.out.vcf : 
//     Channel.empty()

// ch_skip_filtering = params.run_survivor_filter ? 
//     Channel.empty() : 
//     survivor_merge_workflow.out.vcf
//  --run_survivor_filter false

// if close around if params. ......{}
    survivor_filter_workflow(
    survivor_merge_workflow.out.vcf,
    mosdepth_workflow.out.lowcov_bed,
    params.min_sv_size
)

// Get mosdepth directory 
    // ch_spectre_mosdepth = mosdepth_workflow.out.regions_bed
    //     .map { meta, bed_files ->
    //         // Get the parent directory containing all mosdepth output files
    //         def mosdepth_dir = bed_files[0].parent
    //         tuple(meta, mosdepth_dir)
    //     }
// Input parameters - reference genome and tandem repeat file
    ch_spectre_reference = Channel
    .fromPath(params.spectre_snv_vcf, checkIfExists: true)
    .map { vcf_file -> 
        def sample_id = vcf_file.baseName.split('_')[0]
        return [id: sample_id]
    }
    .combine(Channel.fromPath(params.spectre_fasta_file, checkIfExists: true))
    .map { meta, fasta -> tuple(meta, fasta) }

    

// Run SPECTRE CNV calling
   // Run SPECTRE CNV calling
    spectre_cnv_workflow(
        params.spectre_mosdepth,
        ch_spectre_reference,
        params.spectre_snv_vcf,
        params.spectre_metadata,
        params.spectre_blacklist
    )



}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/