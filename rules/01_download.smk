rule download_samples:
    output:
        r1=f"{RAW_DIR}/{{sample}}_1.fastq.gz",
        r2=f"{RAW_DIR}/{{sample}}_2.fastq.gz",
        stats=f"{RAW_DIR}/{{sample}}_data_stats.txt",
    log:
        "logs/download/{sample}.log"
    conda:
        "../envs/01_download.yaml"
    threads:
        config["params"]["download-threads"]
    resources:
        mem_mb=2000
    shell:
        """
        set -euo pipefail

        # Prefetch the SRA object first
        prefetch {wildcards.sample} --output-directory {RAW_DIR}/sra_cache > {log} 2>&1

        # Validate checksums against NCBI's stored values before conversion
        vdb-validate {RAW_DIR}/sra_cache/{wildcards.sample}/{wildcards.sample}.sra >> {log} 2>&1

        # Convert to FASTQ, --split-3 keeps _1/_2 counts consistent
        fasterq-dump {RAW_DIR}/sra_cache/{wildcards.sample}/{wildcards.sample}.sra \
            --split-3 --threads {threads} --outdir {RAW_DIR} >> {log} 2>&1

        # Compress paired reads
        gzip {RAW_DIR}/{wildcards.sample}_1.fastq >> {log} 2>&1
        gzip {RAW_DIR}/{wildcards.sample}_2.fastq >> {log} 2>&1

        # Handle orphan/unpaired reads if --split-3 produced any
        if [ -f {RAW_DIR}/{wildcards.sample}.fastq ]; then
            gzip {RAW_DIR}/{wildcards.sample}.fastq >> {log} 2>&1
        fi

        # Per-sample stats, one header, no shared-file race condition
        seqkit stats {output.r1} {output.r2} > {output.stats}
        """
    

rule download_reference:
    output:
        genome=f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        integrity=f"{REF_DIR}/reference_integrity.txt"
    log:
        "logs/download/reference.log"
    conda:
        "../envs/01_download.yaml"
    threads:
        config["params"]["download-threads"]
    resources:
        mem_mb=2000
    shell:
        """
        set -euo pipefail

        mkdir -p {REF_DIR}/known_sites

        # Download the reference archive and Ensembl's checksum manifest
        wget -O {REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
            https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz >> {log} 2>&1
        wget -O {REF_DIR}/CHECKSUMS_dna \
            https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/CHECKSUMS >> {log} 2>&1

        # Verify the archive against Ensembl's published md5 before using it
        expected=$(grep "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" {REF_DIR}/CHECKSUMS_dna | awk '{{print $1}}')
        actual=$(md5sum {REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz | awk '{{print $1}}')
        if [ "$expected" != "$actual" ]; then
            echo "reference checksum mismatch: expected $expected, got $actual" | tee -a {log}
            exit 1
        fi
        echo "verified md5 $actual against Ensembl CHECKSUMS" > {output.integrity}

        # Decompress to the FASTA consumed downstream
        gunzip -f {REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

        # Download known variant sites (for BQSR)
        wget -O {REF_DIR}/known_sites/Homo_sapiens_assembly38.dbsnp138.vcf \
            https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf >> {log} 2>&1
        wget -O {REF_DIR}/known_sites/Homo_sapiens_assembly38.dbsnp138.vcf.idx \
            https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx >> {log} 2>&1
        wget -O {REF_DIR}/known_sites/Homo_sapiens_assembly38.known_indels.vcf.gz \
            https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz >> {log} 2>&1
        wget -O {REF_DIR}/known_sites/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi \
            https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi >> {log} 2>&1
        """