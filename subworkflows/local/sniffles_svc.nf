#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SNIFFLES } from '../../modules/nf-core/sniffles/main.nf'

workflow sniffles_workflow {
    
    take:
    input         // tuple(val(meta), path(bam), path(bai))
    fasta         // tuple(val(meta), path(fasta))
    tandem_file   // tuple(val(meta), path(bed))
    vcf_output    // val(true)
    snf_output    // val(true)

    main:
    SNIFFLES(
        input,
        fasta,
        tandem_file,
        vcf_output,
        snf_output
    )

    emit:
    vcf     = SNIFFLES.out.vcf
    tbi     = SNIFFLES.out.tbi
    snf     = SNIFFLES.out.snf }
