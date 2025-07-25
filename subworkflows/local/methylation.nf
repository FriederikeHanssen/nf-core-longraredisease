

include { MODKIT_PILEUP } from '../../modules/nf-core/modkit/pileup/main'
include { MODKIT_PILEUP as MODKIT_PILEUP_CPG } from '../../modules/nf-core/modkit/pileup/main'

workflow methyl_subworkflow {

    take:
    bam_bai    // channel: [ val(meta), path(bam), path(bai), path(bed) ]
    fasta       // channel: [ val(meta2), path(fasta) ] - optional
    bed         // channel: [ val(meta3), path(bed) ] - optional

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Run MODKIT PILEUP
    //
    MODKIT_PILEUP (
        bam_bai,
        fasta,
        bed
    )

    MODKIT_PILEUP_CPG (
        bam_bai,
        fasta,
        bed
    )

    ch_versions = ch_versions.mix(MODKIT_PILEUP.out.versions.first())

    emit:
    bed         = MODKIT_PILEUP.out.bed     // channel: [ val(meta), path(bed) ]
    bedgraph    = MODKIT_PILEUP.out.bedgraph // channel: [ val(meta), path(bedgraph) ]
    bed_cpg     = MODKIT_PILEUP_CPG.out.bed     // channel: [ val(meta), path(bed) ]
    bedgraph_cpg = MODKIT_PILEUP_CPG.out.bedgraph // channel: [ val(meta), path(bedgraph) 
    versions = ch_versions               // channel: [ versions.yml ]
}
