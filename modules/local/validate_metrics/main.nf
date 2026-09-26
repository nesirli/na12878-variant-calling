process VALIDATE_METRICS {
    tag "$meta.id"
    label 'process_single'

    container 'docker://quay.io/biocontainers/gawk:5.4.1'

    input:
    tuple val(meta), path(isec_dir)

    output:
    tuple val(meta), path("${meta.id}.validation.txt"), emit: report
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # bcftools isec output: 0000 = calls only, 0001 = truth only, 0002 = both
    YOUR_ONLY=\$(awk '!/^#/ {n++} END {print n+0}' ${isec_dir}/0000.vcf)
    TRUTH_ONLY=\$(awk '!/^#/ {n++} END {print n+0}' ${isec_dir}/0001.vcf)
    CONCORDANT=\$(awk '!/^#/ {n++} END {print n+0}' ${isec_dir}/0002.vcf)

    SENSITIVITY=\$(awk -v c="\$CONCORDANT" -v t="\$TRUTH_ONLY" \\
        'BEGIN { if (c + t > 0) printf "%.3f", c / (c + t); else print "NA" }')
    PRECISION=\$(awk -v c="\$CONCORDANT" -v y="\$YOUR_ONLY" \\
        'BEGIN { if (c + y > 0) printf "%.3f", c / (c + y); else print "NA" }')

    {
        echo "=== Validation: ${meta.id} ==="
        echo "Your calls only:   \$YOUR_ONLY"
        echo "Truth only:        \$TRUTH_ONLY"
        echo "Both (concordant): \$CONCORDANT"
        echo "Sensitivity:       \$SENSITIVITY"
        echo "Precision:         \$PRECISION"
    } | tee ${prefix}.validation.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk --version | head -1 | sed 's/.*Awk //; s/,.*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "stub validation" > ${prefix}.validation.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: 5.4.1
    END_VERSIONS
    """
}
