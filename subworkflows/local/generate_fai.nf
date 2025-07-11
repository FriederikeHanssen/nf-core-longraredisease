include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx/main'


workflow generate_fai_subworkflow {

    take:
    ch_fasta // channel: [meta, fasta] - input FASTA files
    sizes // boolean: whether to generate sizes file

    main:

    // Generate FASTA index (FAI) files
    SAMTOOLS_FAIDX(ch_fasta, sizes)

    ch_fai = SAMTOOLS_FAIDX.out.fai

    emit:
    ch_fai // channel: [meta, fasta.fai] - indexed FASTA files
}