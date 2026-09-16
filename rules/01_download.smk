rule download_samples:
    output:
        r1=f"{RAW_DIR}/{{sample}}_1.fastq.gz",
        r2=f"{RAW_DIR}/{{sample}}_2.fastq.gz",
        stats=f"{RAW_DIR}/data_stats.txt",
    log:
        "log/download/{sample}.log"
    conda:
        "../env/01_download.yaml"
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
        """