// This workflow uses CUTESV to call structural variants (SVs) from BAM files.
// It takes BAM files (with optional FASTA reference) and outputs a VCF file with SV calls.
// It also compresses and indexes the VCF file using bgzip and tabix.

include { CUTESV } from '../../modules/nf-core/cutesv/main.nf'
// include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_CUTESV} from '../../modules/nf-core/tabix/bgziptabix/main.nf'

workflow cutesv_sv_subworkflow {

    take:
    input         // tuple(val(meta), path(bam), path(bai))
    fasta         // tuple(val(meta), path(fasta))

    main:
    // 1. Call SVs with CUTESV
    CUTESV(input, fasta)

    // 2. Compress + index the VCF
    // CUTESV.out.vcf.map { meta, vcf -> tuple(meta, vcf) }
    //               .set { ch_vcf_to_index }

    // TABIX_BGZIPTABIX_CUTESV(ch_vcf_to_index)

    emit:
    vcf         = CUTESV.out.vcf
    // vcf_gz      = TABIX_BGZIPTABIX_CUTESV.out.gz_tbi
    versions    = CUTESV.out.versions
}
