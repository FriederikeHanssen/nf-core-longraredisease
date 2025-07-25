process FILTERCOV_SV {
    tag "$meta.id"
    label 'process_medium'
    
    container "community.wave.seqera.io/library/bcftools_pip_confargparse:4f3c18aa8341a070"
    
    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(mosdepth_summary)
    tuple val(meta3), path(target_bed)
    val chromosome_codes
    val min_read_support
    val min_read_support_limit
    val filter_pass
    
    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: filterbycov_vcf
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_name = "${prefix}"
    def ctgs = chromosome_codes.join(',')
    def ctgs_filter = "--contigs ${ctgs}"
    def pass_filter_arg = filter_pass ? "--filter_pass" : "--no-filter_pass"
    
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

    # Use input VCF directly (already compressed with index)
    input_vcf="${vcf}"
    
    # Verify index exists and is accessible
    if [[ ! -f "${tbi}" ]]; then
        echo "Error: Index file ${tbi} not found"
        exit 1
    fi

    # Extract average depth from mosdepth summary
    AVG_DEPTH=\$(awk '\$1 == "total" {print \$4}' ${mosdepth_summary})
    
    # Generate filtering command with PASS filter option
    get_filter_calls_command.py \\
        --bcftools_threads ${task.cpus} \\
        --target_bedfile \$target_bed_arg \\
        --vcf \$input_vcf \\
        --depth_summary ${mosdepth_summary} \\
        --min_read_support ${min_read_support} \\
        --min_read_support_limit ${min_read_support_limit} \\
        ${ctgs_filter} \\
        ${pass_filter_arg} \\
        ${args} > filter_command.sh
    
    # Execute filtering and compress output with custom naming
    bash filter_command.sh | bcftools view -O z -o ${output_name}.vcf.gz
    
    # Index the output VCF
    tabix -p vcf ${output_name}.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/bcftools //g')
        tabix: \$(tabix --version 2>&1 | head -n1 | sed 's/tabix (htslib) //g')
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "covFiltered"
    def output_name = "${prefix}_${suffix}"
    """
    touch ${output_name}.vcf.gz
    touch ${output_name}.vcf.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/bcftools //g')
        tabix: \$(tabix --version 2>&1 | head -n1 | sed 's/tabix (htslib) //g')
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}