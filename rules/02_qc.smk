rule raw_fastqc:
    input:
        r1=f"{RAW_DIR}/{{sample}}_1.fastq.gz",
        r2=f"{RAW_DIR}/{{sample}}_2.fastq.gz",
    output:
        html1=f"{QC_DIR}/raw/{{sample}}_1_fastqc.html",
        html2=f"{QC_DIR}/raw/{{sample}}_2_fastqc.html",
    log:
        "logs/qc/raw/{sample}.log"
    conda:
        "../envs/02_qc.yaml"
    container:
        "docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    threads:
        config["params"]["qc-threads"]
    resources:
        mem_mb=2000
    shell:
        """
        exec 2> {log}
        set -x
        set -euo pipefail

        mkdir -p {QC_DIR}/raw
        fastqc {input.r1} {input.r2} -o {QC_DIR}/raw --threads {threads}
        """

rule raw_multi_qc:
    input:
        f"{QC_DIR}/raw/{sample}_1_fastqc.html",
        f"{QC_DIR}/raw/{sample}_2_fastqc.html",
    output:
        f"{QC_DIR}/raw/multiqc_report.html"
    log:
        "logs/qc/raw/multiqc.log"
    conda:
        "../envs/02_qc.yaml"
    container:
        "docker://quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1"
    threads:
        1
    resources:
        mem_mb=2000
    shell:
        """
        exec 2> {log}
        set -x
        set -euo pipefail

        multiqc {QC_DIR}/raw -o {QC_DIR}/raw
        """