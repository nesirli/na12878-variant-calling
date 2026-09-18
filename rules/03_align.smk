rule bwa_index:
    input:
        f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    output:
        multiext(f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
                 ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "logs/align/bwa_index.log"
    conda:
        "../envs/03_align.yaml"
    container:
        "docker://quay.io/biocontainers/bwa:0.7.19--h577a1d6_1"
    threads:
        1
    resources:
        mem_mb=16000
    shell:
        """
        bwa index {input} > {log} 2>&1
        """


rule bwa_mem:
    input:
        r1=f"{RAW_DIR}/{{sample}}_1.fastq.gz",
        r2=f"{RAW_DIR}/{{sample}}_2.fastq.gz",
        ref=f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        ref_idx=multiext(f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
                         ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        pipe(f"{ALIGN_DIR}/{{sample}}.unsorted.sam")
    log:
        "logs/align/{sample}.bwa.log"
    conda:
        "../envs/03_align.yaml"
    container:
        "docker://quay.io/biocontainers/bwa:0.7.19--h577a1d6_1"
    threads:
        config["params"]["bwa-threads"]
    resources:
        mem_mb=12000
    shell:
        """
        exec 2> {log}
        set -x
        set -euo pipefail

        bwa mem -t {threads} \
            -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA\\tLB:lib1' \
            {input.ref} \
            {input.r1} {input.r2} > {output}
        """


rule samtools_sort:
    input:
        sam=f"{ALIGN_DIR}/{{sample}}.unsorted.sam"
    output:
        bam=f"{ALIGN_DIR}/{{sample}}.sorted.bam",
        bai=f"{ALIGN_DIR}/{{sample}}.sorted.bam.bai",
        stats=f"{ALIGN_DIR}/{{sample}}_samtools_stats.txt"
    log:
        "logs/align/{sample}.samtools.log"
    conda:
        "../envs/03_align.yaml"
    container:
        "docker://quay.io/biocontainers/samtools:1.24--h9dcdb79_1"
    threads:
        config["params"]["sort-threads"]
    params:
        sort_mem_mb=config["params"]["sort-mem-mb"],
    resources:
        mem_mb=6000
    shell:
        """
        exec 2> {log}
        set -x
        set -euo pipefail

        samtools sort -@ {threads} -m {params.sort_mem_mb}M -o {output.bam} {input.sam}

        samtools index {output.bam}

        samtools flagstat {output.bam} > {output.stats}
        """