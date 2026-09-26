process ANN_FILTER {
    tag "$meta.id"
    label 'process_single'

    container 'docker://quay.io/biocontainers/gawk:5.4.1'

    input:
    tuple val(meta), path(vcf)
    path(filter_awk)

    output:
    tuple val(meta), path("${meta.id}.high_impact.vcf"), emit: high_impact
    tuple val(meta), path("${meta.id}.missense.vcf")   , emit: missense
    tuple val(meta), path("${meta.id}.annotation_summary.txt"), emit: summary
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Most frequent variant effects (ANN field, 2nd sub-field)
    grep -v '^#' ${vcf} \\
        | sed 's/.*ANN=//' \\
        | cut -d'|' -f2 \\
        | sort | uniq -c | sort -rn | head -15 \\
        > ${prefix}.annotation_summary.txt || true

    # High-impact variants (frameshift, stop gained, splice)
    awk -v mode=high -f ${filter_awk} ${vcf} > ${prefix}.high_impact.vcf
    # Missense variants
    awk -v mode=missense -f ${filter_awk} ${vcf} > ${prefix}.missense.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk --version | head -1 | sed 's/.*Awk //; s/,.*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.high_impact.vcf ${prefix}.missense.vcf ${prefix}.annotation_summary.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: 5.4.1
    END_VERSIONS
    """
}
