
include { SURVIVOR_MERGE } from '../../modules/nf-core/survivor/merge/main.nf'
include { SURVIVOR_EXTRACT_INTERVALS } from '../../modules/local/extract_survivor_bed/main.nf'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_SNIFFLES} from '../../modules/nf-core/bedtools/intersect/main.nf'

workflow multi_caller_sv_filter_subworkflow {

    take:
    survivor_vcfs        // Channel: tuple(meta, List[VCF file]) 
    sv_vcf_sniffles      // Channel: tuple(meta, path(sv_vcf_sniffles))
    chrom_sizes          // Channel: tuple(meta, path(chrom_sizes)) - optional
    max_dist             // int
    min_callers          // int
    type_flag            // int
    strand_flag          // int
    est_dist             // int
    min_size             // int

    main:

    survivor_vcfs
        .map { meta, files ->
            tuple(meta, files)
        }
        .set { ch_vcf_input }


    SURVIVOR_MERGE(
        ch_vcf_input,
        max_dist,
        min_callers,
        type_flag,
        strand_flag,
        est_dist,
        min_size
    )

    SURVIVOR_EXTRACT_INTERVALS(
        SURVIVOR_MERGE.out.vcf
    )


    // Prepare channels for filtering VCF
    
    ch_input_bedtools_sniffles = sv_vcf_sniffles
    .combine(SURVIVOR_EXTRACT_INTERVALS.out.bed)
    .filter { vcf_meta, vcf_file, bed_meta, bed_file -> 
        vcf_meta.id == bed_meta.id 
    }
    .map { vcf_meta, vcf_file, bed_meta, bed_file ->
        tuple(vcf_meta, vcf_file, bed_file)
    }
    
    // Filter each VCF using bedtools intersect
    BEDTOOLS_INTERSECT_SNIFFLES(ch_input_bedtools_sniffles, chrom_sizes)

    emit:
    survivor_vcf          = SURVIVOR_MERGE.out.vcf
    bed                   = SURVIVOR_EXTRACT_INTERVALS.out.bed
    filtered_vcf_sniffles = BEDTOOLS_INTERSECT_SNIFFLES.out.intersect
}