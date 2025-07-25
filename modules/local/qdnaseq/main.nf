process QDNASEQ {
    tag "$meta.id"
    label 'process_medium'
    
    container "docker.io/nourmahfel1/runqdnaseq:0.0.4"

    input:
    tuple val(meta), path(bam), path(bai)
    val genome_build
    val bin_size
    val method
    val cutoff
    val cutoff_del
    val cutoff_loss
    val cutoff_gain
    val cellularity

    output:
    tuple val(meta), path("*_calls.vcf")        , emit: calls_vcf
    tuple val(meta), path("*_segs.vcf")         , emit: segs_vcf
    tuple val(meta), path("*_raw_calls.vcf")    , emit: raw_calls_vcf, optional: true
    tuple val(meta), path("*_raw_segs.vcf")     , emit: raw_segs_vcf, optional: true
    tuple val(meta), path("*.pdf")              , emit: plots, optional: true
    tuple val(meta), path("*.png")              , emit: images, optional: true
    tuple val(meta), path("*.txt")              , emit: tables, optional: true
    tuple val(meta), path("*.bed")              , emit: bed, optional: true
    tuple val(meta), path("*.seg")              , emit: seg, optional: true
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = genome_build ?: "hg38"
    def binsize = bin_size ?: 500
    def analysis_method = method ?: "cutoff"
    def cn_cutoff = cutoff ?: 0.5
    def cn_cutoff_del = cutoff_del ?: 0.5
    def cn_cutoff_loss = cutoff_loss ?: 1.5
    def cn_cutoff_gain = cutoff_gain ?: 2.5
    def cn_cellularity = cellularity ?: 1
    
    """
    # Run QDNAseq analysis
    Rscript /app/run_qdnaseq.r \\
        --bam ${bam} \\
        --out_prefix ${prefix} \\
        --reference ${reference} \\
        --binsize ${binsize} \\
        --method ${analysis_method} \\
        --cutoff ${cn_cutoff} \\
        --cutoffDEL ${cn_cutoff_del} \\
        --cutoffLOSS ${cn_cutoff_loss} \\
        --cutoffGAIN ${cn_cutoff_gain} \\
        --cellularity ${cn_cellularity} \\
        ${args}
    
    # Check if raw VCF files were created and rename them
    if [ -f "${prefix}_calls.vcf" ]; then
        mv ${prefix}_calls.vcf raw_calls.vcf
    fi
    
    if [ -f "${prefix}_segs.vcf" ]; then
        mv ${prefix}_segs.vcf raw_segs.vcf
    fi
    
    # Fix VCF files using the provided Python script
    if [ -f "raw_calls.vcf" ]; then
        python /app/fix_1491_vcf.py -i raw_calls.vcf -o ${prefix}_calls.vcf --sample_id ${prefix}
        mv raw_calls.vcf ${prefix}_raw_calls.vcf
    fi
    
    if [ -f "raw_segs.vcf" ]; then
        python /app/fix_1491_vcf.py -i raw_segs.vcf -o ${prefix}_segs.vcf --sample_id ${prefix}
        mv raw_segs.vcf ${prefix}_raw_segs.vcf
    fi
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        runqdnaseq: "0.0.2"
        r-base: "4.2.2"
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_calls.vcf
    touch ${prefix}_segs.vcf
    touch ${prefix}_raw_calls.vcf
    touch ${prefix}_raw_segs.vcf
    touch ${prefix}_plots.pdf
    touch ${prefix}_analysis.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        runqdnaseq: "0.0.1"
        qdnaseq: "1.38.0"
        r-base: "4.2.2"
        container: "docker.io/nourmahfel1/runqdnaseq:0.0.1"
    END_VERSIONS
    """
}