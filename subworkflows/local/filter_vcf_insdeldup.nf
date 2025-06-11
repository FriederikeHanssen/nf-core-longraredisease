// Optionally, this subworkflow can be used to filter a VCF file for insertions, deletions, and duplications.
// process filterBenchmarkVcf {
//     label "wf_human_sv"
//     cpus 2
//     memory 4.GB
//     input:
//         tuple val(xam_meta), path(calls_vcf)
//     output:
//         tuple val(xam_meta), path("benchmarkCalls.vcf.gz"), path("benchmarkCalls.vcf.gz.tbi")
//     script:
//     """
//     zcat $calls_vcf \
//     | bcftools view -i '(SVTYPE = \"INS\" || SVTYPE = \"DEL\" || SVTYPE = \"DUP\")' \
//     | bgziptabix benchmarkCalls.vcf.gz
//     """
// }
