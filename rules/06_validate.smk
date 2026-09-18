rule download_truth_set:
    output:
        vcf=f"{config['truth']['dir']}/{config['truth']['filename']}",
        tbi=f"{config['truth']['dir']}/{config['truth']['filename']}.tbi",
    log:
        "logs/validate/download_truth.log"
    conda:
        "../envs/01_download.yaml"
    container:
        "docker://quay.io/biocontainers/wget:1.25.0"
    params:
        url=config["truth"]["url"],
        dir=config["truth"]["dir"],
        filename=config["truth"]["filename"],
    shell:
        """
        exec 2> {log}
        set -x
        set -euo pipefail

        mkdir -p {params.dir}
        wget -c -O {output.vcf} {params.url}/{params.filename}
        wget -c -O {output.tbi} {params.url}/{params.filename}.tbi
        """


rule prepare_truth_set:
    input:
        vcf=f"{config['truth']['dir']}/{config['truth']['filename']}",
        tbi=f"{config['truth']['dir']}/{config['truth']['filename']}.tbi",
    output:
        vcf=f"{VALIDATION_DIR}/truth.ensembl.vcf.gz",
        tbi=f"{VALIDATION_DIR}/truth.ensembl.vcf.gz.tbi",
    log:
        "logs/validate/prepare_truth.log"
    conda:
        "../envs/06_validate.yaml"
    container:
        "docker://quay.io/biocontainers/bcftools:1.24--h118bc1c_2"
    shell:
        """
        exec 2> {log}
        set -x
        set -euo pipefail

        # GIAB uses UCSC-style names (chr20); the Ensembl reference uses 20.
        # Rename contigs so bcftools isec can intersect the two call sets.
        # This truth set covers chromosomes 1-22 only.
        for i in {{1..22}}; do
            printf 'chr%s\t%s\n' "$i" "$i"
        done > {VALIDATION_DIR}/chr_rename.txt

        bcftools annotate --rename-chrs {VALIDATION_DIR}/chr_rename.txt \
            {input.vcf} -Oz -o {output.vcf}
        bcftools index -t {output.vcf}
        """


rule validate:
    input:
        calls=f"{ANALYSIS_DIR}/{{sample}}.final.vcf.gz",
        calls_index=f"{ANALYSIS_DIR}/{{sample}}.final.vcf.gz.tbi",
        truth=f"{VALIDATION_DIR}/truth.ensembl.vcf.gz",
        truth_index=f"{VALIDATION_DIR}/truth.ensembl.vcf.gz.tbi",
    output:
        report=f"{VALIDATION_DIR}/{{sample}}.validation.txt",
    log:
        "logs/validate/{sample}.validate.log"
    params:
        interval=config["params"]["variant-interval"],
    conda:
        "../envs/06_validate.yaml"
    container:
        "docker://quay.io/biocontainers/bcftools:1.24--h118bc1c_2"
    shell:
        """
        exec 2> {log}
        set -x
        set -euo pipefail

        isec_dir={VALIDATION_DIR}/{wildcards.sample}.isec
        rm -rf "$isec_dir"

        # Restrict to the called interval (empty = whole genome)
        region=""
        if [ -n "{params.interval}" ]; then
            region="-r {params.interval}"
        fi

        # 0000 = calls only, 0001 = truth only, 0002 = concordant sites
        bcftools isec -p "$isec_dir" $region {input.calls} {input.truth}

        YOUR_ONLY=$(awk '!/^#/ {{n++}} END {{print n+0}}' "$isec_dir/0000.vcf")
        TRUTH_ONLY=$(awk '!/^#/ {{n++}} END {{print n+0}}' "$isec_dir/0001.vcf")
        CONCORDANT=$(awk '!/^#/ {{n++}} END {{print n+0}}' "$isec_dir/0002.vcf")

        SENSITIVITY=$(awk -v c="$CONCORDANT" -v t="$TRUTH_ONLY" \
            'BEGIN {{ if (c + t > 0) printf "%.3f", c / (c + t); else print "NA" }}')
        PRECISION=$(awk -v c="$CONCORDANT" -v y="$YOUR_ONLY" \
            'BEGIN {{ if (c + y > 0) printf "%.3f", c / (c + y); else print "NA" }}')

        {{
            echo "=== Validation: {wildcards.sample} ==="
            echo "Your calls only:   $YOUR_ONLY"
            echo "Truth only:        $TRUTH_ONLY"
            echo "Both (concordant): $CONCORDANT"
            echo "Sensitivity:       $SENSITIVITY"
            echo "Precision:         $PRECISION"
        }} | tee {output.report}
        """
