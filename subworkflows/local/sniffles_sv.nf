#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SNIFFLES } from '../../modules/nf-core/sniffles/main.nf'
include { GUNZIP as GUNZIP_SNIFFLES } from '../../modules/nf-core/gunzip/main.nf'

workflow sniffles_sv_subworkflow {
    
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

    // Only run GUNZIP if VCF files were produced
    ch_vcf_uncompressed = Channel.empty()
    if (vcf_output) {
        GUNZIP_SNIFFLES(SNIFFLES.out.vcf)
        ch_vcf_uncompressed = GUNZIP_SNIFFLES.out.gunzip
    }

    emit:
    vcf         = ch_vcf_uncompressed     // Uncompressed VCF files
    // vcf_gz      = SNIFFLES.out.vcf        // Original compressed VCF files
    // tbi         = SNIFFLES.out.tbi        // Tabix index files
    // snf         = SNIFFLES.out.snf        // SNF files
    versions    = SNIFFLES.out.versions   // Version information
}