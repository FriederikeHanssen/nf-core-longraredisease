
include { NANOPLOT } from '../../modules/nf-core/nanoplot/main'


workflow nanoplot_subworkflow {

    take:
    ch_bam        // Channel: tuple(meta, path(bam), path(bai))

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Run NANOPLOT
    //
    
    NANOPLOT (
        ch_bam
    )

    ch_versions = ch_versions.mix(NANOPLOT.out.versions.first())

    emit:
    html_report = NANOPLOT.out.html  // Channel: tuple(meta, path(html_report))
    versions = ch_versions                   // Channel: [ versions.yml ]
}