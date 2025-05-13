/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { sniffles_workflow } from '../subworkflows/local/sniffles_svc.nf'
include { cutesv_workflow    } from '../subworkflows/local/cutesv_svc.nf'
include { svim_workflow      } from '../subworkflows/local/svim_svc.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow nanoraredx {

    ch_input_bam = Channel
        .fromPath(params.bam_file, checkIfExists: true)
        .map { bam ->
            def bai = file("${bam}.bai")
            def meta = [ id: bam.getBaseName().replaceFirst(/\.bam$/, '') ]
            tuple(meta, bam, bai)
        }

    ch_fasta = Channel
        .fromPath(params.fasta_file, checkIfExists: true)
        .map { fasta -> tuple([id: "ref"], fasta) }

    ch_trf = Channel
        .fromPath(params.tandem_file, checkIfExists: true)
        .map { bed -> tuple([id: "trf"], bed) }

    // Parallel execution
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