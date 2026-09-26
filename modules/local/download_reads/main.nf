process DOWNLOAD_READS {
    tag "$meta.id"
    label 'process_single'

    container 'docker://quay.io/biocontainers/wget:1.25.0'

    input:
    val(meta)

    output:
    tuple val(meta), path("${meta.id}_*.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.fastq.gz")  , emit: orphan, optional: true
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def base = task.ext.base_url ?: "${params.ena_url}"
    def id   = meta.id
    def md5_r1     = params.ena_md5_r1     ? "echo '${params.ena_md5_r1}  ${id}_1.fastq.gz' | md5sum -c -"     : "true"
    def md5_r2     = params.ena_md5_r2     ? "echo '${params.ena_md5_r2}  ${id}_2.fastq.gz' | md5sum -c -"     : "true"
    def md5_orphan = params.ena_md5_orphan ? "echo '${params.ena_md5_orphan}  ${id}.fastq.gz' | md5sum -c -"   : "true"
    """
    # ENA occasionally drops long transfers and briefly answers 403/5xx: retry + resume
    WGET="wget -c $args --tries=100 --waitretry=60 --timeout=120 --read-timeout=120 \\
        --retry-connrefused --retry-on-host-error --no-http-keep-alive \\
        --retry-on-http-error=403,408,429,500,502,503,504"

    \$WGET -O ${id}_1.fastq.gz '${base}/${id}_1.fastq.gz'
    \$WGET -O ${id}_2.fastq.gz '${base}/${id}_2.fastq.gz'
    \$WGET -O ${id}.fastq.gz   '${base}/${id}.fastq.gz'

    # Verify against ENA's published checksums when provided
    ${md5_r1}
    ${md5_r2}
    ${md5_orphan}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo \$(wget --version 2>&1) | sed 's/^.*GNU Wget //; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_1.fastq.gz ${meta.id}_2.fastq.gz ${meta.id}.fastq.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: 1.25.0
    END_VERSIONS
    """
}
