process CLAIR3 {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "docker.io/hkubal/clair3"

    input:
    tuple val(meta), path(bam), path(bai), val(packaged_model), path(user_model), val(platform)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(index)

    output:
    tuple val(meta), path("*merge_output.vcf.gz"),            emit: vcf
    tuple val(meta), path("*merge_output.vcf.gz.tbi"),        emit: tbi
    tuple val(meta), path("*full_alignment.vcf.gz"),          emit: full_vcf
    tuple val(meta), path("*full_alignment.vcf.gz.tbi"),      emit: full_tbi
    tuple val(meta), path("*pileup.vcf.gz"),                  emit: pileup_vcf
    tuple val(meta), path("*pileup.vcf.gz.tbi"),              emit: pileup_tbi
    tuple val(meta), path("*phased_merge_output.vcf.gz"),     emit: phased_vcf, optional: true
    tuple val(meta), path("*phased_merge_output.vcf.gz.tbi"), emit: phased_tbi, optional: true
    path "versions.yml",                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def model = ""
    if (!user_model) {
        if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
            model = "\${CONDA_PREFIX}/opt/models/${packaged_model}"
        }
        else {
            model = "/opt/models/$packaged_model"
        }
    }
    if (!packaged_model) {
        model = "$user_model"
    }
    if (packaged_model && user_model) {
        error "Two models specified $user_model and $packaged_model, specify one of them."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    run_clair3.sh \\
        --bam_fn=$bam \\
        --ref_fn=$reference \\
        --threads=$task.cpus \\
        --sample_name=${meta.id} \\
        --output=. \\
        --platform=$platform \\
        --model=$model \\
        $args
    
    mv merge_output.vcf.gz ${prefix}_merge_output.vcf.gz
    mv merge_output.vcf.gz.tbi ${prefix}_merge_output.vcf.gz.tbi

    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(run_clair3.sh  --version |& sed '1!d ; s/Clair3 v//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.phased_merge_output.vcf.gz
    touch ${prefix}.phased_merge_output.vcf.gz.tbi
    echo "" | gzip > ${prefix}.merge_output.vcf.gz
    touch ${prefix}.merge_output.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(run_clair3.sh --version |& sed '1!d ; s/Clair3 v//')
    END_VERSIONS
    """
}
