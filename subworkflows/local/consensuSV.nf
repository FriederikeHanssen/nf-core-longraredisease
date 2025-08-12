include { SURVIVOR_MERGE                         } from '../../modules/nf-core/survivor/merge/main.nf'
include { SURVIVOR_FIX                            } from '../../modules/local/consensuSV/survivor_fixVCF/main.nf'
include { EXTRACT_BED                       } from '../../modules/local/consensuSV/extract_bed/main.nf'
include { BCFTOOLS_MERGE  as BCFTOOLS_MERGE_SV    } from '../../modules/nf-core/bcftools/merge/main.nf'
include { TABIX_BGZIPTABIX as BGZIP_CONSENSUSV      } from '../../modules/nf-core/tabix/bgziptabix/main.nf'
include { TRUVARI_COLLAPSE                        } from '../../modules/local/truvari/collapse/main.nf'
include { MERGE_SURV_TRUVARI                 } from '../../modules/local/consensuSV/merge_surv_truvari/main.nf'


workflow consensuSV_subworkflow {
    take:
    survivor_vcfs        // Channel: tuple(meta, List[VCF file]) 
    bcftools_vcfs       // Channel: tuple(meta, List[VCF.gz file], List[TBI file])

    main:
    ch_versions = Channel.empty()

    // Step 1: Run SURVIVOR_MERGE
    SURVIVOR_MERGE(
        survivor_vcfs,
        params.max_distance_breakpoints ?: 1000,
        params.min_supporting_callers ?: 2,
        params.account_for_type ?: true,
        params.account_for_sv_strands ?: true,
        params.estimate_distanced_by_sv_size ?: false,
        params.min_sv_size ?: 30
    )
    SURVIVOR_FIX(SURVIVOR_MERGE.out.vcf) 
    EXTRACT_BED(SURVIVOR_FIX.out.vcf)

    ch_versions = ch_versions.mix(SURVIVOR_MERGE.out.versions)

    // Step 3: bcftools merge
    BCFTOOLS_MERGE_SV(
        bcftools_vcfs,
        [[:], []],
        [[:], []],
        [[:], []]
    )

    BGZIP_CONSENSUSV(BCFTOOLS_MERGE_SV.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE_SV.out.versions)

    // TRUVARI_COLLAPSE
    TRUVARI_COLLAPSE(
        BGZIP_CONSENSUSV.out.gz_tbi, 
        params.refdist, 
        params.pctsim, 
        params.pctseq,
        params.passonly, 
        params.dup_to_ins
    )
    ch_versions = ch_versions.mix(TRUVARI_COLLAPSE.out.versions)



    MERGE_SURV_TRUVARI(
        TRUVARI_COLLAPSE.out.merged_vcf,
        EXTRACT_BED.out.bed
    )

    emit:
    vcf = MERGE_SURV_TRUVARI.out.vcf
    versions = ch_versions
}