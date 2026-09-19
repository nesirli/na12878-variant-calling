import os

# SnpEff resolves a relative -dataDir against its own install directory, so we
# always pass an absolute path based on the workflow directory.
SNPEFF_DATA_DIR = os.path.join(
    workflow.basedir, config["snpeff"]["data-dir"]
)
SNPEFF_DB_DIR = os.path.join(SNPEFF_DATA_DIR, config["snpeff"]["database"])
SNPEFF_MARKER = os.path.join(SNPEFF_DB_DIR, "snpEffectPredictor.bin")


rule snpeff_download:
    output:
        db=SNPEFF_MARKER,
    log:
        "logs/annotate/snpeff_download.log"
    conda:
        "../envs/05_annotate.yaml"
    container:
        "docker://quay.io/biocontainers/snpeff:5.1--hdfd78af_4"
    params:
        database=config["snpeff"]["database"],
        data_dir=SNPEFF_DATA_DIR,
    shell:
        """
        exec 2> {log}
        set -x
        set -euo pipefail

        mkdir -p {params.data_dir}
        snpEff download -v -dataDir {params.data_dir} {params.database}
        """


rule annotate:
    input:
        vcf=f"{ANALYSIS_DIR}/{{sample}}.final.vcf.gz",
        db=SNPEFF_MARKER,
    output:
        annotated=f"{ANNOTATION_DIR}/{{sample}}.annotated.vcf",
        annotation_summary=f"{ANNOTATION_DIR}/{{sample}}.annotation_summary.txt",
        high_impact=f"{ANNOTATION_DIR}/{{sample}}.high_impact.vcf",
        missense=f"{ANNOTATION_DIR}/{{sample}}.missense.vcf",
    log:
        "logs/annotate/{sample}.annotate.log"
    conda:
        "../envs/05_annotate.yaml"
    container:
        "docker://quay.io/biocontainers/snpeff:5.1--hdfd78af_4"
    params:
        database=config["snpeff"]["database"],
        data_dir=SNPEFF_DATA_DIR,
    shell:
        """
        exec 2> {log}
        set -x
        set -euo pipefail

        # Annotate with SnpEff (it reads bgzipped VCF input directly)
        # -Xmx is parsed by the snpEff wrapper and forwarded to the JVM; loading
        # the GRCh38.105 database needs more than the 1g default heap.
        snpEff -Xmx16g ann -noStats -dataDir {params.data_dir} {params.database} \
            {input.vcf} > {output.annotated}

        # Summary of the most frequent variant effects (ANN field, 2nd sub-field)
        grep -v '^#' {output.annotated} \
            | sed 's/.*ANN=//' \
            | cut -d'|' -f2 \
            | sort | uniq -c | sort -rn | head -15 \
            > {output.annotation_summary} || true

        # Extract high-impact variants (frameshift, stop gained, splice).
        # The bioconda snpeff package ships snpEff but not SnpSift, so parse
        # the ANN field directly (see scripts/filter_ann.awk).
        awk -v mode=high -f scripts/filter_ann.awk {output.annotated} \
            > {output.high_impact}
        echo "High-impact variants: $(grep -vc '^#' {output.high_impact})"

        # Extract missense variants
        awk -v mode=missense -f scripts/filter_ann.awk {output.annotated} \
            > {output.missense}
        echo "Missense variants: $(grep -vc '^#' {output.missense})"
        """
