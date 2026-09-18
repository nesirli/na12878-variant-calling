# Load the configuration file
configfile: "config/config.yaml"

# Access the variables
RAW_DIR = config["directories"]["raw"]
REF_DIR = config["directories"]["reference"]
QC_DIR = config["directories"]["qc"]
ALIGN_DIR = config["directories"]["align"]

SAMPLES = ['ERR250949']

rule all:
    input:
        f"{QC_DIR}/raw/multiqc_report.html",
        expand(f"{ALIGN_DIR}/{{sample}}.sorted.bam", sample=SAMPLES),
        expand(f"{ALIGN_DIR}/{{sample}}.sorted.bam.bai", sample=SAMPLES),
        expand(f"{ALIGN_DIR}/{{sample}}_samtools_stats.txt", sample=SAMPLES),
        expand(f"{RAW_DIR}/{{sample}}_data_stats.txt", sample=SAMPLES)

include: "rules/01_download.smk"
include: "rules/02_qc.smk"
include: "rules/03_align.smk"
#include: "rules/04_analysis.smk"
# "rules/05_variants.smk"
# "rules/06_annotate.smk"
# "rules/07_validate.smk"