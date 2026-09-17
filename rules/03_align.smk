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
    threads:
        config["params"]["align-threads"]
    resources:
        mem_mb=16000
    shell:
        """
        bwa index {input} > {log} 2>&1
        """


rule align:
    input:
        r1=f"{RAW_DIR}/{{sample}}_1.fastq.gz",
        r2=f"{RAW_DIR}/{{sample}}_2.fastq.gz",
        ref=f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        ref_idx=multiext(f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
                         ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        bam=f"{ALIGN_DIR}/{{sample}}.sorted.bam",
        bai=f"{ALIGN_DIR}/{{sample}}.sorted.bam.bai",
        stats=f"{ALIGN_DIR}/{{sample}}_samtools_stats.txt"
    log:
        "logs/align/{sample}.log"
    conda:
        "../envs/03_align.yaml"
    threads:
        config["params"]["align-threads"]
    resources:
        mem_mb=16000
    shell:
        """
        exec 2> {log}
        set -x
        set -euo pipefail

        bwa mem -t {threads} \
            -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA\\tLB:lib1' \
            {input.ref} \
            {input.r1} {input.r2} | \
            samtools sort -@ 4 -o {output.bam} -

        samtools index {output.bam}

        samtools flagstat {output.bam} > {output.stats}

        """