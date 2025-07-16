// This workflow uses CUTESV to call structural variants (SVs) from BAM files.
// It takes BAM files (with optional FASTA reference) and outputs a VCF file with SV calls.
// It also compresses and indexes the VCF file using bgzip and tabix.

include { CUTESV } from '../../modules/nf-core/cutesv/main.nf'


workflow cutesv_sv_subworkflow {

    take:
    input         // tuple(val(meta), path(bam), path(bai))
    fasta         // tuple(val(meta), path(fasta))

    main:

    CUTESV(input, fasta)



    emit:
    vcf         = CUTESV.out.vcf
    versions    = CUTESV.out.versions
}
