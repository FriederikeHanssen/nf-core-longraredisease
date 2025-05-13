#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SVIM } from '../../modules/local/svim/main.nf'

workflow svim_workflow {

    take:
    input         // tuple(val(meta), path(bam), path(bai))
    fasta         // tuple(val(meta), path(fasta))

    main:
    SVIM(
        input,
        fasta
    )
    // Emit the output files
    emit:
    vcf     = SVIM.out.vcf
    versions = SVIM.out.versions
}

