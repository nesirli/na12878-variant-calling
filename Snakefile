# Load the configuration file
configfile: "config/config.yaml"

# Access the variables
RAW_DIR = config["directories"]["raw"]
REF_DIR = config["directories"]["reference"]
QC_DIR = config["directories"]["qc"]
ALIGN_DIR = config["directories"]["align"]
ANALYSIS_DIR = config["directories"]["analysis"]
ANNOTATION_DIR = config["directories"]["annotation"]
VALIDATION_DIR = config["directories"]["validation"]

sample = config["sample"]["id"]

rule all:
    input:
        f"{QC_DIR}/raw/multiqc_report.html",
        f"{ALIGN_DIR}/{sample}.sorted.bam",
        f"{ALIGN_DIR}/{sample}.sorted.bam.bai",
        f"{ALIGN_DIR}/{sample}_samtools_stats.txt",
        f"{RAW_DIR}/{sample}_data_stats.txt",
        f"{ANALYSIS_DIR}/{sample}.final.vcf.gz",
        f"{ANNOTATION_DIR}/{sample}.annotated.vcf",
        f"{VALIDATION_DIR}/{sample}.validation.txt"

include: "rules/01_download.smk"
include: "rules/02_qc.smk"
include: "rules/03_align.smk"
include: "rules/04_analysis.smk"
include: "rules/05_annotate.smk"
include: "rules/06_validate.smk"