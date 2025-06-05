// Subworkflow for filtering SVs using SURVIVOR with optional low coverage BED mask

include { SURVIVOR_FILTER } from '../../modules/nf-core/survivor/filter/main.nf'

workflow survivor_filter_workflow {
    take:
    survivor_merged_vcf  // channel: tuple(val(meta), path(vcf))
    lowcov_bed            // channel: tuple(val(meta), path(bed))
    minsv                 // val: e.g. 50
    maxsv                 // val: e.g. -1 (no max)
    minallelefreq         // val: e.g. 0.0
    minnumreads           // val: e.g. -1 (disable read support filter)

    main:
    // Combine inputs into one tuple per sample
    ch_filter_input = survivor_merged_vcf
        .join(lowcov_bed, by: 0, remainder: true)
        .map { meta, vcf, bed ->
            def bed_file = bed ?: []  // Use empty list if no BED file
            tuple(meta, vcf, bed_file)
        }


    SURVIVOR_FILTER(
        ch_filter_input,
        minsv,
        maxsv,
        minallelefreq,
        minnumreads
    )

    emit:
    vcf      = SURVIVOR_FILTER.out.vcf
    versions = SURVIVOR_FILTER.out.versions
}
