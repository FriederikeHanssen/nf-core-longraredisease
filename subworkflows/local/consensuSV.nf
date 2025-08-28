include { SURVIVOR_MERGE                         } from '../../modules/nf-core/survivor/merge/main.nf'
include { SURVIVOR_VCFTOBED                      } from '../../modules/local/survivor_vcftobed/main.nf'
include { BCFTOOLS_MERGE  as BCFTOOLS_MERGE_SV   } from '../../modules/nf-core/bcftools/merge/main.nf'
include { TABIX_BGZIPTABIX as BGZIP_CONSENSUSV   } from '../../modules/nf-core/tabix/bgziptabix/main.nf'
include { TRUVARI_COLLAPSE                       } from '../../modules/local/truvari/collapse/main.nf'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_SV      } from '../../modules/nf-core/bcftools/sort/main.nf'
include { TABIX_BGZIPTABIX as TRUVARI_GZ         } from '../../modules/nf-core/tabix/bgziptabix/main.nf'

workflow consensuSV_subworkflow {
    take:
    survivor_vcfs        // Channel: tuple(meta, List[VCF file]) 
    bcftools_vcfs       // Channel: tuple(meta, List[VCF.gz file], List[TBI file])
    use_survivor_bed     // Boolean: true to use SURVIVOR, false to skip it

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
    

    if (use_survivor_bed) {
        SURVIVOR_MERGE(
            survivor_vcfs,
            params.max_distance_breakpoints,
            params.min_supporting_callers,
            params.account_for_type,
            params.account_for_sv_strands,
            params.estimate_distanced_by_sv_size,
            params.min_sv_size
        )

       
        SURVIVOR_VCFTOBED(SURVIVOR_MERGE.out.vcf)
        ch_vcf = SURVIVOR_MERGE.out.vcf
        ch_bed = SURVIVOR_VCFTOBED.out.bed
        ch_versions = ch_versions.mix(SURVIVOR_MERGE.out.versions)
    } else {
        
        ch_vcf = Channel.empty()
        ch_bed = Channel.empty()
    }
    
    ch_vcf_normalized = BGZIP_CONSENSUSV.out.gz_tbi
        .map { meta, vcf, tbi ->
            def normalized_meta = [id: meta.id]  
            tuple(normalized_meta, vcf, tbi)
        }
    

    if (use_survivor_bed) {
        ch_bed_normalized = ch_bed
            .map { meta, bed ->
                def normalized_meta = [id: meta.id]  
                tuple(normalized_meta, bed)
            }
        
        // Join the normalized channels
        ch_combined_input = ch_vcf_normalized
            .join(ch_bed_normalized, by: 0)  
            .map { meta, vcf, tbi, bed ->      //this is redundant
                tuple(meta, vcf, tbi, bed)
            }
    } else {
        // Create dummy BED channel for TRUVARI_COLLAPSE
        ch_combined_input = ch_vcf_normalized
            .map { meta, vcf, tbi ->
                def dummy_bed = file("NO_FILE")
                tuple(meta, vcf, tbi, dummy_bed)
            }
    }

    TRUVARI_COLLAPSE(
        ch_combined_input,               // tuple val(meta), path(vcf), path(tbi), path(bed)
        params.refdist,          // val(refdist)
        params.pctsim,            // val(pctsim)
        params.pctseq,            // val(pctseq)
        params.passonly,         // val(passonly)
        params.dup_to_ins     // val(dup_to_ins)
    )

    ch_versions = ch_versions.mix(TRUVARI_COLLAPSE.out.versions)

    // Step 6: Compress and index TRUVARI output
    BCFTOOLS_SORT_SV(TRUVARI_COLLAPSE.out.merged_vcf)
    TRUVARI_GZ(BCFTOOLS_SORT_SV.out.vcf)

    emit:
    vcf = TRUVARI_COLLAPSE.out.merged_vcf
    vcf_gz = TRUVARI_GZ.out.gz_tbi
    collapsed_vcf = TRUVARI_COLLAPSE.out.collapsed_vcf
    survivor_vcf = ch_vcf
    survivor_bed = ch_bed
    versions = ch_versions
}