include { SURVIVOR_MERGE                         } from '../../modules/nf-core/survivor/merge/main.nf'
include { SURVIVOR_VCFTOBED                      } from '../../modules/local/survivor_vcftobed/main.nf'
include { BCFTOOLS_MERGE  as BCFTOOLS_MERGE_SV   } from '../../modules/nf-core/bcftools/merge/main.nf'
include { TABIX_BGZIPTABIX as BGZIP_CONSENSUSV   } from '../../modules/nf-core/tabix/bgziptabix/main.nf'
include { TRUVARI_COLLAPSE                       } from '../../modules/local/truvari/collapse/main.nf'
include { TABIX_BGZIPTABIX as TRUVARI_GZ         } from '../../modules/nf-core/tabix/bgziptabix/main.nf'

workflow consensuSV_subworkflow {
    take:
    survivor_vcfs        // Channel: tuple(meta, List[VCF file]) 
    bcftools_vcfs       // Channel: tuple(meta, List[VCF.gz file], List[TBI file])

    main:
    ch_versions = Channel.empty()

    BCFTOOLS_MERGE_SV(
        bcftools_vcfs,
        [[:], []],
        [[:], []],
        [[:], []]
    )

    BGZIP_CONSENSUSV(BCFTOOLS_MERGE_SV.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE_SV.out.versions)
    
    SURVIVOR_MERGE(
        survivor_vcfs,
        params.max_distance_breakpoints ?: 1000,
        params.min_supporting_callers ?: 2,
        params.account_for_type ?: true,
        params.account_for_sv_strands ?: true,
        params.estimate_distanced_by_sv_size ?: false,
        params.min_sv_size ?: 30
    )

    // Step 2: Extract BED from SURVIVOR output
    SURVIVOR_VCFTOBED(SURVIVOR_MERGE.out.vcf)
    ch_vcf = SURVIVOR_MERGE.out.vcf
    ch_bed = SURVIVOR_VCFTOBED.out.bed
    ch_versions = ch_versions.mix(SURVIVOR_MERGE.out.versions)
    
    // DEBUG: Check the meta objects before joining
    BGZIP_CONSENSUSV.out.gz_tbi.view { "BGZIP meta: $it" }
    ch_bed.view { "BED meta: $it" }
    
    // **FIX: Normalize meta objects before joining**
    ch_vcf_normalized = BGZIP_CONSENSUSV.out.gz_tbi
        .map { meta, vcf, tbi ->
            def normalized_meta = [id: meta.id]  // Keep only the ID
            tuple(normalized_meta, vcf, tbi)
        }
    
    ch_bed_normalized = ch_bed
        .map { meta, bed ->
            def normalized_meta = [id: meta.id]  // Keep only the ID
            tuple(normalized_meta, bed)
        }
    
    // DEBUG: Check normalized channels
    ch_vcf_normalized.view { "VCF normalized: $it" }
    ch_bed_normalized.view { "BED normalized: $it" }
    
    // Join the normalized channels
    ch_combined_input = ch_vcf_normalized
        .join(ch_bed_normalized, by: 0)  // Join by meta (index 0)
        .map { meta, vcf, tbi, bed ->
            tuple(meta, vcf, tbi, bed)
        }
    
    // DEBUG: Check final combined channel
    ch_combined_input.view { "Combined input: $it" }

    TRUVARI_COLLAPSE(
        ch_combined_input,               // tuple val(meta), path(vcf), path(tbi), path(bed)
        params.refdist ?: 1000,          // val(refdist)
        params.pctsim ?: 0.7,            // val(pctsim)
        params.pctseq ?: 0.7,            // val(pctseq)
        params.passonly ?: true,         // val(passonly)
        params.dup_to_ins ?: false       // val(dup_to_ins)
    )

    ch_versions = ch_versions.mix(TRUVARI_COLLAPSE.out.versions)

    // Step 6: Compress and index TRUVARI output
    TRUVARI_GZ(TRUVARI_COLLAPSE.out.merged_vcf)

    emit:
    vcf = TRUVARI_COLLAPSE.out.merged_vcf
    vcf_gz = TRUVARI_GZ.out.gz_tbi
    collapsed_vcf = TRUVARI_COLLAPSE.out.collapsed_vcf
    survivor_vcf = ch_vcf
    survivor_bed = ch_bed
    versions = ch_versions
}