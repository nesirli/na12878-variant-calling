process DOWNLOAD_REFERENCE {
    tag 'GRCh38'
    label 'process_single'

    container 'docker://quay.io/biocontainers/wget:1.25.0'

    input:
    val(meta)

    output:
    tuple val(meta), path("${params.reference_name}.fa"), emit: fasta
    tuple val(meta), path("known_sites/Homo_sapiens_assembly38.dbsnp138.vcf"), path("known_sites/Homo_sapiens_assembly38.dbsnp138.vcf.idx"), emit: dbsnp
    tuple val(meta), path("known_sites/Homo_sapiens_assembly38.known_indels.vcf.gz"), path("known_sites/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"), emit: indels
    path "reference_integrity.txt", emit: integrity
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p known_sites

    # Reference FASTA + Ensembl checksum manifest (resumable)
    wget -c $args -O ${params.reference_name}.fa.gz '${params.fasta_url}'
    wget -O CHECKSUMS_dna '${params.checksums_url}'

    # Verify against the Ensembl BSD-sum manifest before decompressing
    expected=\$(awk '\$3 == "${params.reference_name}.fa.gz" {print \$1, \$2}' CHECKSUMS_dna)
    actual=\$(sum ${params.reference_name}.fa.gz | awk '{print \$1, \$2}')
    if [ "\$expected" != "\$actual" ]; then
        echo "reference checksum mismatch: expected \$expected, got \$actual" >&2
        exit 1
    fi
    echo "verified sum \$actual against Ensembl CHECKSUMS" > reference_integrity.txt

    gunzip -f ${params.reference_name}.fa.gz

    # Known sites for BQSR (dbsnp + known indels)
    wget -c -O known_sites/Homo_sapiens_assembly38.dbsnp138.vcf '${params.dbsnp_url}'
    wget -c -O known_sites/Homo_sapiens_assembly38.dbsnp138.vcf.idx '${params.dbsnp_idx_url}'
    wget -c -O known_sites/Homo_sapiens_assembly38.known_indels.vcf.gz '${params.indels_url}'
    wget -c -O known_sites/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi '${params.indels_idx_url}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo \$(wget --version 2>&1) | sed 's/^.*GNU Wget //; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    touch ${params.reference_name}.fa
    mkdir -p known_sites
    touch known_sites/Homo_sapiens_assembly38.dbsnp138.vcf
    touch known_sites/Homo_sapiens_assembly38.dbsnp138.vcf.idx
    touch known_sites/Homo_sapiens_assembly38.known_indels.vcf.gz
    touch known_sites/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
    echo "stub" > reference_integrity.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: 1.25.0
    END_VERSIONS
    """
}
