#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CUTESV } from '../../modules/nf-core/cutesv/main.nf'

workflow cutesv_workflow {

    take:
    input
    fasta

    main:
    CUTESV(
        input,
        fasta)

    emit:
    vcf = CUTESV.out.vcf
}