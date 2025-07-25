include { SURVIVOR_MERGE                          } from '../../modules/nf-core/survivor/merge/main.nf'
include { SURVIVOR_FIX                            } from '../../modules/local/survivor_fix_vcf/main.nf'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_SV       } from '../../modules/nf-core/bcftools/sort/main.nf'
include {TABIX_BGZIPTABIX as BGZIP_SV             } from '../../modules/nf-core/tabix/bgziptabix/main.nf'
// include { BCFTOOLS_MERGE as BCFTOOLS_MERGE_SV     }   from '../../modules/nf-core/bcftools/merge/main.nf'
// include { TRUVARI_COLLAPSE                        } from '../../modules/local/truvari/collapse/main.nf'
// include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_SV } from '../../modules/nf-core/bedtools/intersect/main.nf'
// include { BCFTOOLS_NORM                              } from '../../modules/nf-core/bcftools/norm/main'    


// change the naming to bgziptabix sv
workflow consensuSV_subworkflow {

    take:
    survivor_vcfs        // Channel: tuple(meta, List[VCF file]) 
    max_dist             // int
    min_callers          // int
    type_flag            // int
    strand_flag          // int
    est_dist             // int
    min_size             // int

    main:

    // Ensure the input channel is properly structured for SURVIVOR_MERGE
    ch_vcf_input = survivor_vcfs
        .map { meta, files ->
            // Make sure files is a list and meta is properly structured
            tuple(meta, files)
        }

    SURVIVOR_MERGE(
        ch_vcf_input,
        max_dist,
        min_callers,
        type_flag,
        strand_flag,
        est_dist,
        min_size
    )

    
    SURVIVOR_FIX(SURVIVOR_MERGE.out.vcf) 

    BCFTOOLS_SORT_SV(SURVIVOR_FIX.out.vcf)
    BGZIP_SV(BCFTOOLS_SORT_SV.out.vcf)

    ch_vcf_gz = BGZIP_SV.out.gz_tbi
            .map { meta, vcf_gz, tbi -> 
                tuple(meta, vcf_gz) 
            }

        

    emit:
    vcf             = BCFTOOLS_SORT_SV.out.vcf // channel: [meta, vcf]
    vcf_gz          = ch_vcf_gz // channel: [meta, vcf.gz]
    versions        = SURVIVOR_MERGE.out.versions

}