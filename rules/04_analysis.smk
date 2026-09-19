rule samtools_faidx:
    input:
        ref=f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    output:
        fai=f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"
    log:
        "logs/analysis/samtools_faidx.log"
    conda:
        "../envs/04_analysis.yaml"
    container:
        "docker://quay.io/biocontainers/samtools:1.24--h9dcdb79_1"
    shell:
        """
        exec 2> {log}
        set -x
        set -euo pipefail

        samtools faidx {input.ref}
        """


rule gatk_create_dict:
    input:
        ref=f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    output:
        dict=f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.dict"
    log:
        "logs/analysis/gatk_create_dict.log"
    conda:
        "../envs/04_analysis.yaml"
    container:
        "docker://quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0"
    shell:
        """
        exec 2> {log}
        set -x
        set -euo pipefail

        gatk CreateSequenceDictionary -R {input.ref} -O {output.dict}
        """


rule mark_duplicates:
    input:
        sorted_bam=f"{ALIGN_DIR}/{{sample}}.sorted.bam",
    output:
        metrics=f"{ANALYSIS_DIR}/{{sample}}.dedup_metrics.txt",
        dedup_bam=f"{ANALYSIS_DIR}/{{sample}}.dedup.bam",
        dedup_bai=f"{ANALYSIS_DIR}/{{sample}}.dedup.bai",
    log:
        "logs/analysis/{sample}.mark_duplicates.log"
    conda:
        "../envs/04_analysis.yaml"
    container:
        "docker://quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0"
    shell:
        """
        exec 2> {log}
        set -x
        set -euo pipefail

        gatk MarkDuplicates \
            -I {input.sorted_bam} \
            -O {output.dedup_bam} \
            -M {output.metrics} \
            --REMOVE_DUPLICATES false \
            --CREATE_INDEX true
        """


rule prepare_known_sites:
    input:
        dbsnp=f"{REF_DIR}/known_sites/Homo_sapiens_assembly38.dbsnp138.vcf",
        indels=f"{REF_DIR}/known_sites/Homo_sapiens_assembly38.known_indels.vcf.gz",
    output:
        dbsnp=f"{REF_DIR}/known_sites/Homo_sapiens_assembly38.dbsnp138.ensembl.vcf.gz",
        dbsnp_idx=f"{REF_DIR}/known_sites/Homo_sapiens_assembly38.dbsnp138.ensembl.vcf.gz.tbi",
        indels=f"{REF_DIR}/known_sites/Homo_sapiens_assembly38.known_indels.ensembl.vcf.gz",
        indels_idx=f"{REF_DIR}/known_sites/Homo_sapiens_assembly38.known_indels.ensembl.vcf.gz.tbi",
    log:
        "logs/analysis/prepare_known_sites.log"
    conda:
        "../envs/04_analysis.yaml"
    container:
        "docker://quay.io/biocontainers/bcftools:1.24--h118bc1c_2"
    shell:
        """
        exec 2> {log}
        set -x
        set -euo pipefail

        # The Broad known-sites VCFs use UCSC-style contig names (chr20); the
        # Ensembl reference uses 20. Rename so GATK sees overlapping contigs.
        map={REF_DIR}/known_sites/chr_rename.txt
        for i in {{1..22}}; do
            printf 'chr%s\t%s\n' "$i" "$i"
        done > "$map"
        printf 'chrX\tX\nchrY\tY\nchrM\tMT\n' >> "$map"

        bcftools annotate --threads 8 --rename-chrs "$map" {input.dbsnp} -Oz -o {output.dbsnp}
        bcftools index -t {output.dbsnp}
        bcftools annotate --threads 8 --rename-chrs "$map" {input.indels} -Oz -o {output.indels}
        bcftools index -t {output.indels}
        """


rule build_recalibration_model:
    input:
        dedup_bam=f"{ANALYSIS_DIR}/{{sample}}.dedup.bam",
        ref=f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        dict=f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.dict",
        fai=f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai",
        dbsnp=f"{REF_DIR}/known_sites/Homo_sapiens_assembly38.dbsnp138.ensembl.vcf.gz",
        dbsnp_idx=f"{REF_DIR}/known_sites/Homo_sapiens_assembly38.dbsnp138.ensembl.vcf.gz.tbi",
        indels=f"{REF_DIR}/known_sites/Homo_sapiens_assembly38.known_indels.ensembl.vcf.gz",
        indels_idx=f"{REF_DIR}/known_sites/Homo_sapiens_assembly38.known_indels.ensembl.vcf.gz.tbi",
    output:
        recal_table=f"{ANALYSIS_DIR}/{{sample}}.recal.table",
    log:
        "logs/analysis/{sample}.build_recalibration_model.log"
    conda:
        "../envs/04_analysis.yaml"
    container:
        "docker://quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0"
    shell:
        """
        exec 2> {log}
        set -x
        set -euo pipefail

        gatk BaseRecalibrator \
            -I {input.dedup_bam} \
            -R {input.ref} \
            --known-sites {input.dbsnp} \
            --known-sites {input.indels} \
            -O {output.recal_table}
        """


rule apply_recalibration:
    input:
        dedup_bam=f"{ANALYSIS_DIR}/{{sample}}.dedup.bam",
        recal_table=f"{ANALYSIS_DIR}/{{sample}}.recal.table",
        ref=f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
    output:
        recal_bam=f"{ANALYSIS_DIR}/{{sample}}.recal.bam",
        recal_bai=f"{ANALYSIS_DIR}/{{sample}}.recal.bai",
    log:
        "logs/analysis/{sample}.apply_recalibration.log"
    conda:
        "../envs/04_analysis.yaml"
    container:
        "docker://quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0"
    shell:
        """
        exec 2> {log}
        set -x
        set -euo pipefail

        gatk ApplyBQSR \
            -I {input.dedup_bam} \
            -R {input.ref} \
            --bqsr-recal-file {input.recal_table} \
            -O {output.recal_bam} \
            --create-output-bam-index true
        """


rule call_variants:
    input:
        recal_bam=f"{ANALYSIS_DIR}/{{sample}}.recal.bam",
        ref=f"{REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
    output:
        raw_vcf=f"{ANALYSIS_DIR}/{{sample}}.raw.vcf.gz",
    log:
        "logs/analysis/{sample}.call_variants.log"
    params:
        interval=config["params"]["variant-interval"],
    conda:
        "../envs/04_analysis.yaml"
    container:
        "docker://quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0"
    shell:
        """
        exec 2> {log}
        set -x
        set -euo pipefail

        # Restrict to an interval when one is configured (empty = whole genome)
        interval=""
        if [ -n "{params.interval}" ]; then
            interval="-L {params.interval}"
        fi

        gatk HaplotypeCaller \
            -I {input.recal_bam} \
            -R {input.ref} \
            $interval \
            -O {output.raw_vcf}
        """


rule filter_snps:
    input:
        raw_vcf=f"{ANALYSIS_DIR}/{{sample}}.raw.vcf.gz",
    output:
        filtered_snps=f"{ANALYSIS_DIR}/{{sample}}.filtered_snps.vcf.gz",
    log:
        "logs/analysis/{sample}.filter_snps.log"
    conda:
        "../envs/04_analysis.yaml"
    container:
        "docker://quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0"
    shell:
        """
        exec 2> {log}
        set -x
        set -euo pipefail

        gatk SelectVariants \
            -V {input.raw_vcf} \
            --select-type-to-include SNP \
            -O {ANALYSIS_DIR}/{wildcards.sample}.raw_snps.vcf.gz

        gatk VariantFiltration -V {ANALYSIS_DIR}/{wildcards.sample}.raw_snps.vcf.gz \
            --filter-expression "QD < 2.0" --filter-name "LowQD" \
            --filter-expression "MQ < 40.0" --filter-name "LowMQ" \
            --filter-expression "FS > 60.0" --filter-name "HighFS" \
            --filter-expression "SOR > 3.0" --filter-name "HighSOR" \
            -O {output.filtered_snps}
        """


rule filter_indels:
    input:
        raw_vcf=f"{ANALYSIS_DIR}/{{sample}}.raw.vcf.gz",
    output:
        filtered_indels=f"{ANALYSIS_DIR}/{{sample}}.filtered_indels.vcf.gz",
    log:
        "logs/analysis/{sample}.filter_indels.log"
    conda:
        "../envs/04_analysis.yaml"
    container:
        "docker://quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0"
    shell:
        """
        exec 2> {log}
        set -x
        set -euo pipefail

        gatk SelectVariants \
            -V {input.raw_vcf} \
            --select-type-to-include INDEL \
            -O {ANALYSIS_DIR}/{wildcards.sample}.raw_indels.vcf.gz

        gatk VariantFiltration -V {ANALYSIS_DIR}/{wildcards.sample}.raw_indels.vcf.gz \
            --filter-expression "QD < 2.0" --filter-name "LowQD" \
            --filter-expression "FS > 200.0" --filter-name "HighFS" \
            -O {output.filtered_indels}
        """


rule combine_variants:
    input:
        filtered_snps=f"{ANALYSIS_DIR}/{{sample}}.filtered_snps.vcf.gz",
        filtered_indels=f"{ANALYSIS_DIR}/{{sample}}.filtered_indels.vcf.gz",
    output:
        final_vcf=f"{ANALYSIS_DIR}/{{sample}}.final.vcf.gz",
        final_vcf_index=f"{ANALYSIS_DIR}/{{sample}}.final.vcf.gz.tbi",
        snp_count=f"{ANALYSIS_DIR}/{{sample}}.snps.txt",
        indel_count=f"{ANALYSIS_DIR}/{{sample}}.indels.txt",
        ts_tv=f"{ANALYSIS_DIR}/{{sample}}.ts_tv.txt",
    log:
        "logs/analysis/{sample}.combine_variants.log"
    conda:
        "../envs/04_analysis.yaml"
    container:
        "docker://quay.io/biocontainers/bcftools:1.24--h118bc1c_2"
    shell:
        """
        exec 2> {log}
        set -x
        set -euo pipefail

        # Keep only PASS variants
        bcftools view -f PASS -Oz \
            {input.filtered_snps} \
            -o {ANALYSIS_DIR}/{wildcards.sample}.pass_snps.vcf.gz
        bcftools view -f PASS -Oz \
            {input.filtered_indels} \
            -o {ANALYSIS_DIR}/{wildcards.sample}.pass_indels.vcf.gz

        # Index the per-type VCFs; `concat -a` requires indexed inputs
        bcftools index -t {ANALYSIS_DIR}/{wildcards.sample}.pass_snps.vcf.gz
        bcftools index -t {ANALYSIS_DIR}/{wildcards.sample}.pass_indels.vcf.gz

        # Merge SNPs and indels (overlapping positions are allowed; output is sorted)
        bcftools concat -a -Oz \
            -o {output.final_vcf} \
            {ANALYSIS_DIR}/{wildcards.sample}.pass_snps.vcf.gz \
            {ANALYSIS_DIR}/{wildcards.sample}.pass_indels.vcf.gz
        bcftools index -t {output.final_vcf}

        # Summary
        bcftools view -v snps {output.final_vcf} | grep -vc '^#' > {output.snp_count} || true
        bcftools view -v indels {output.final_vcf} | grep -vc '^#' > {output.indel_count} || true
        bcftools stats {output.final_vcf} | grep '^TSTV' > {output.ts_tv} || true
        """
