process FILTER_SV_CALLS {
    tag "$meta.id"
    label 'process_medium'
    
    container "community.wave.seqera.io/library/bcftools_pip_confargparse:4f3c18aa8341a070"
    
    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(mosdepth_summary)
    tuple val(meta3), path(target_bed)
    val chromosome_codes
    val min_read_support
    val min_read_support_limit
    
    output:
    tuple val(meta), path("*.filtered.vcf"), emit: filtered_vcf
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ctgs = chromosome_codes.join(',')
    def ctgs_filter = "--contigs ${ctgs}"
    """
    # Filter bed file for callable regions - handle gzipped files directly
    if [[ "${target_bed}" != "OPTIONAL_FILE" ]]; then
        if [[ "${target_bed}" == *.gz ]]; then
            zcat ${target_bed} | awk '\$4 == "CALLABLE" || \$4 == "HIGH_COVERAGE"' > callable_regions.bed
        else
            awk '\$4 == "CALLABLE" || \$4 == "HIGH_COVERAGE"' ${target_bed} > callable_regions.bed
        fi
        target_bed_arg="callable_regions.bed"
    else
        target_bed_arg=""
    fi

    # Compress and index input VCF if needed
    if [[ "${vcf}" != *.gz ]]; then
        bcftools view -O z ${vcf} > input.vcf.gz
        tabix -p vcf input.vcf.gz
        input_vcf="input.vcf.gz"
    else
        input_vcf="${vcf}"
        if [[ ! -f "${vcf}.tbi" ]]; then
            tabix -p vcf ${vcf}
        fi
    fi

    # Extract average depth from mosdepth summary
    AVG_DEPTH=\$(awk '\$1 == "total" {print \$4}' ${mosdepth_summary})
    
    # Generate filtering command
    get_filter_calls_command.py \\
        --bcftools_threads ${task.cpus} \\
        --target_bedfile \$target_bed_arg \\
        --vcf \$input_vcf \\
        --depth_summary ${mosdepth_summary} \\
        --min_read_support ${min_read_support} \\
        --min_read_support_limit ${min_read_support_limit} \\
        ${ctgs_filter} \\
        ${args} > filter_command.sh
    
    # Execute filtering
    bash filter_command.sh > ${prefix}.filtered.vcf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/bcftools //g')
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.filtered.vcf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/bcftools //g')
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}