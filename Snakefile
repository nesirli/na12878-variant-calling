# Load the configuration file
configfile: "config/config.yaml"

# Access the variables
RAW_DIR = config["directories"]["raw"]
REF_DIR = config["directories"]["reference"]
QC_DIR = config["directories"]["qc"]

SAMPLES = ['ERR250949']

rule all:
    input:
        expand(f"{RAW_DIR}/{{sample}}_1.fastq.gz", sample=SAMPLES),
        expand(f"{RAW_DIR}/{{sample}}_2.fastq.gz", sample=SAMPLES),
        expand(f"{RAW_DIR}/{{sample}}_data_stats.txt", sample=SAMPLES),
        f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        f"{REF_DIR}/reference_integrity.txt"

include: "rules/01_download.smk"
# "rules/02_qc.smk"
# "rules/03_align.smk"
# "rules/04_analysis.smk"
# "rules/05_variants.smk"
# "rules/06_annotate.smk"
# "rules/07_validate.smk"