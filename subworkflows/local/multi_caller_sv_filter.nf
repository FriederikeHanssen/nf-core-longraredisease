#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SURVIVOR_MERGE } from '../../modules/nf-core/survivor/merge/main.nf'
include { SURVIVOR_EXTRACT_INTERVALS } from '../../modules/local/extract_survivor_bed/main.nf'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_CUTESV } from '../../modules/nf-core/bedtools/intersect/main.nf'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_SNIFFLES} from '../../modules/nf-core/bedtools/intersect/main.nf'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_SVIM } from '../../modules/nf-core/bedtools/intersect/main.nf'
include { TABIX_BGZIPTABIX as BGZIP_CUTESV } from '../../modules/nf-core/tabix/bgziptabix/main.nf'
include { TABIX_BGZIPTABIX as BGZIP_SNIFFLES } from '../../modules/nf-core/tabix/bgziptabix/main.nf'
include { TABIX_BGZIPTABIX as BGZIP_SVIM } from '../../modules/nf-core/tabix/bgziptabix/main.nf'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_SV} from '../../modules/nf-core/bcftools/concat/main.nf'

workflow multi_caller_sv_filter_subworkflow {

    take:
    survivor_vcfs        // Channel: tuple(meta, List[VCF file])
    sv_vcf_cutesv        // Channel: tuple(meta, path(sv_vcf_cutesv))
    sv_vcf_sniffles      // Channel: tuple(meta, path(sv_vcf_sniffles))
    sv_vcf_svim          // Channel: tuple(meta, path(sv_vcf_svim))
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

    // Prepare channels for filtering each VCF
    ch_input_bedtools_cutesv = sv_vcf_cutesv
            .join(SURVIVOR_EXTRACT_INTERVALS.out.bed, by: 0)
            .map { meta, vcf, bed ->
                def id_parts = meta.id.split('_')
                def sample_id = id_parts[0]
                def updated_meta = meta + [id: "${sample_id}_cutesv"]
                tuple(updated_meta, vcf, bed)
            }
    
    ch_input_bedtools_sniffles = sv_vcf_sniffles
            .join(SURVIVOR_EXTRACT_INTERVALS.out.bed, by: 0)
            .map { meta, vcf, bed ->
                def id_parts = meta.id.split('_')
                def sample_id = id_parts[0]
                def updated_meta = meta + [id: "${sample_id}_sniffles"]
                tuple(updated_meta, vcf, bed)
            }
    
    ch_input_bedtools_svim = sv_vcf_svim
            .join(SURVIVOR_EXTRACT_INTERVALS.out.bed, by: 0)
            .map { meta, vcf, bed ->
                def id_parts = meta.id.split('_')
                def sample_id = id_parts[0]
                def updated_meta = meta + [id: "${sample_id}_svim"]
                tuple(updated_meta, vcf, bed)
            }
    
    // Filter each VCF using bedtools intersect
    BEDTOOLS_INTERSECT_CUTESV(ch_input_bedtools_cutesv, chrom_sizes)
    BEDTOOLS_INTERSECT_SNIFFLES(ch_input_bedtools_sniffles, chrom_sizes)
    BEDTOOLS_INTERSECT_SVIM(ch_input_bedtools_svim, chrom_sizes)

    // Compress and index the filtered VCFs
    BGZIP_CUTESV(BEDTOOLS_INTERSECT_CUTESV.out.intersect)
    BGZIP_SNIFFLES(BEDTOOLS_INTERSECT_SNIFFLES.out.intersect)
    BGZIP_SVIM(BEDTOOLS_INTERSECT_SVIM.out.intersect)

    // Prepare compressed VCFs for concatenation
    ch_filtered_for_concat = BGZIP_CUTESV.out.gz_tbi
        .mix(BGZIP_SNIFFLES.out.gz_tbi)
        .mix(BGZIP_SVIM.out.gz_tbi)
        .map { meta, vcf_gz, tbi ->
            // Extract the base sample ID
            def base_id = meta.id.replaceAll(/_(cutesv|sniffles|svim)$/, '')
            def base_meta = meta + [id: base_id]
            tuple(base_meta, vcf_gz, tbi)
        }
        .groupTuple(by: 0)
        .map { meta, vcf_list, tbi_list ->
            // Sort VCFs and TBIs to ensure consistent order
            def sorted_vcfs = vcf_list.sort()
            def sorted_tbis = tbi_list.sort()
            tuple(meta, sorted_vcfs, sorted_tbis)
        }

    // Concatenate the compressed VCFs
    BCFTOOLS_CONCAT_SV(ch_filtered_for_concat)

    emit:
    survivor_vcf         = SURVIVOR_MERGE.out.vcf
    bed                  = SURVIVOR_EXTRACT_INTERVALS.out.bed
    filtered_vcf_cutesv  = BEDTOOLS_INTERSECT_CUTESV.out.intersect
    filtered_vcf_sniffles = BEDTOOLS_INTERSECT_SNIFFLES.out.intersect
    filtered_vcf_svim    = BEDTOOLS_INTERSECT_SVIM.out.intersect
    concatenated_vcf     = BCFTOOLS_CONCAT_SV.out.vcf
    versions             = SURVIVOR_MERGE.out.versions.mix(BCFTOOLS_CONCAT_SV.out.versions)
}