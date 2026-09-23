# ENA occasionally drops long transfers and briefly answers with 403/5xx.
# Retry aggressively and always resume (-c) so multi-GB downloads survive.
WGET = (
    "wget -c --tries=100 --waitretry=60 --timeout=120 --read-timeout=120 "
    "--retry-connrefused --retry-on-host-error --no-http-keep-alive "
    "--retry-on-http-error=403,408,429,500,502,503,504"
)


rule download_samples:
    output:
        r1=f"{RAW_DIR}/{{sample}}_1.fastq.gz",
        r2=f"{RAW_DIR}/{{sample}}_2.fastq.gz",
        orphan=f"{RAW_DIR}/{{sample}}.fastq.gz",
    log:
        "logs/download/{sample}.log"
    params:
        url=config["sample"]["ena-url"],
        md5_r1=config["sample"]["ena-md5"]["r1"],
        md5_r2=config["sample"]["ena-md5"]["r2"],
        md5_orphan=config["sample"]["ena-md5"]["orphan"],
    conda:
        "../envs/01_download.yaml"
    container:
        "docker://quay.io/biocontainers/wget:1.25.0"
    resources:
        mem_mb=2000
    shell:
        """
        exec 2> {log}
        set -x
        set -euo pipefail

        mkdir -p {RAW_DIR}

        # The ENA mirror serves this run as gzipped FASTQ
        {WGET} -O {output.r1} {params.url}/{wildcards.sample}_1.fastq.gz
        {WGET} -O {output.r2} {params.url}/{wildcards.sample}_2.fastq.gz
        {WGET} -O {output.orphan} {params.url}/{wildcards.sample}.fastq.gz

        # Verify the downloads against ENA's published checksums
        echo "{params.md5_r1}  {output.r1}" | md5sum -c -
        echo "{params.md5_r2}  {output.r2}" | md5sum -c -
        echo "{params.md5_orphan}  {output.orphan}" | md5sum -c -
        """


rule sample_stats:
    input:
        r1=f"{RAW_DIR}/{{sample}}_1.fastq.gz",
        r2=f"{RAW_DIR}/{{sample}}_2.fastq.gz",
    output:
        f"{RAW_DIR}/{{sample}}_data_stats.txt",
    log:
        "logs/download/{sample}_stats.log"
    conda:
        "../envs/01_download.yaml"
    container:
        "docker://quay.io/biocontainers/seqkit:2.13.0--he881be0_0"
    resources:
        mem_mb=2000
    shell:
        """
        exec 2> {log}
        set -x
        set -euo pipefail

        seqkit stats {input.r1} {input.r2} > {output}
        """
    

rule download_reference:
    output:
        genome=f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        integrity=f"{REF_DIR}/reference_integrity.txt",
        dbsnp=f"{REF_DIR}/known_sites/Homo_sapiens_assembly38.dbsnp138.vcf",
        dbsnp_idx=f"{REF_DIR}/known_sites/Homo_sapiens_assembly38.dbsnp138.vcf.idx",
        indels=f"{REF_DIR}/known_sites/Homo_sapiens_assembly38.known_indels.vcf.gz",
        indels_idx=f"{REF_DIR}/known_sites/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"
    log:
        "logs/download/reference.log"
    conda:
        "../envs/01_download.yaml"
    container:
        "docker://quay.io/biocontainers/wget:1.25.0"
    threads:
        config["params"]["download-threads"]
    resources:
        mem_mb=2000
    shell:
        """
        exec 2> {log}
        set -x
        set -euo pipefail

        mkdir -p {REF_DIR}/known_sites

        # Download the reference archive and Ensembl's checksum manifest
        wget -c -O {REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
            https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
        wget -O {REF_DIR}/CHECKSUMS_dna \
            https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/CHECKSUMS

        # Verify the archive against Ensembl's published sum (checksum + block count) before using it
        expected=$(awk '$3 == "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" {{print $1, $2}}' {REF_DIR}/CHECKSUMS_dna)
        actual=$(sum {REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz | awk '{{print $1, $2}}')
        if [ "$expected" != "$actual" ]; then
            echo "reference checksum mismatch: expected $expected, got $actual" | tee -a {log}
            exit 1
        fi
        echo "verified sum $actual against Ensembl CHECKSUMS" > {output.integrity}

        # Decompress to the FASTA consumed downstream
        gunzip -f {REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

        # Download known variant sites (for BQSR); -c resumes partial downloads
        wget -c -O {output.dbsnp} \
            https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
        wget -c -O {output.dbsnp_idx} \
            https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
        wget -c -O {output.indels} \
            https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
        wget -c -O {output.indels_idx} \
            https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
        """