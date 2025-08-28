process SURVIVOR_VCFTOBED {
    tag "$meta.id"
    label 'process_low'
    
    conda "bioconda::bcftools=1.17 bioconda::gawk=5.1.0"
    container "community.wave.seqera.io/library/bcftools_gawk:c2387a9d5226f9b9"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    # Create the awk script file
    cat > process_coords.awk << 'EOF'
function emit(ch,s,e){
    if (s>e){t=s;s=e;e=t}
    s=s-1; 
    if (s<0) s=0;
    if (s<e) printf("%s\\t%d\\t%d\\n", ch, s, e)
}
{
    gsub(/^[ \\t]+/, "", \$0)
    gsub(/[ \\t]+\$/, "", \$0)
    if (\$0=="") next;
    if (match(\$0,/^(chr[0-9A-Za-z_.]+)[_:]([0-9]+)-(chr[0-9A-Za-z_.]+)[_:]([0-9]+)\$/, m)) {
        if (m[1]==m[3]) emit(m[1], m[2], m[4]); 
        next
    }
    if (match(\$0,/^(chr[0-9A-Za-z_.]+)[_:]([0-9]+)-([0-9]+)\$/, m)) {
        emit(m[1], m[2], m[3]); 
        next
    }
}
EOF

    # Run the pipeline in one go
    bcftools query -f '[%CHROM\\t%POS\\t%SAMPLE\\t%CO\\n]' ${vcf} | \\
    awk -F'\\t' 'NF==4 && \$4!="." && \$4!="NA" && \$4!="NAN" && \$4!="0" {print \$4}' | \\
    tr ',' '\\n' | \\
    gawk -f process_coords.awk | \\
    sort -k1,1 -k2,2n > ${prefix}.bed

    # Clean up the awk script
    rm -f process_coords.awk

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        gawk: \$(gawk --version 2>&1 | head -n1 | sed 's/^.*gawk //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        gawk: \$(gawk --version 2>&1 | head -n1 | sed 's/^.*gawk //; s/ .*\$//')
    END_VERSIONS
    """
}