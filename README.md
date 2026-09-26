# na12878-variant-calling

An end-to-end reproducible workflow for calling, filtering, annotating, and validating
clinical genomic variants for the GIAB reference sample **NA12878 / HG001** (ENA run
`ERR194147`), against GRCh38.

This repository contains **two independent implementations** of the same analysis DAG:

| Implementation | Entry point | Status | Docs |
|---|---|---|---|
| Snakemake (original) | `Snakefile`, `rules/`, `config/config.yaml` | on `main` | [docs/README.snakemake.md](docs/README.snakemake.md) |
| Nextflow + nf-core modules | `main.nf`, `workflows/`, `conf/` | on `feat/nf-core-nextflow` | [docs/NEXTFLOW_TUTORIAL.md](docs/NEXTFLOW_TUTORIAL.md) |

The Snakemake workflow remains the reference implementation. The Nextflow port reuses
[curated nf-core modules](https://github.com/nf-core/modules) so each step uses a
community-maintained, containerised process definition, and is designed to run on a
Slurm cluster launched from [Seqera Cloud](https://cloud.seqera.io).

## Analysis stages (common to both)

1. Download reads (ENA) and GRCh38 reference + known sites
2. Read QC (FastQC / MultiQC)
3. Alignment (BWA-MEM) and sorting/indexing (samtools)
4. GATK4 MarkDuplicates, BaseRecalibrator, ApplyBQSR, HaplotypeCaller
5. Variant filtering (GATK) and merging (bcftools)
6. Functional annotation (SnpEff)
7. Concordance against the GIAB v4.2.1 truth set (bcftools isec)

## Snakemake usage (main)

```bash
snakemake --cores 12 --sdm conda --resources mem_mb=45000 --scheduler ilp
```

See [docs/README.snakemake.md](docs/README.snakemake.md) for full instructions.

## Nextflow usage (this branch)

```bash
nextflow run . -profile apptainer,slurm --input assets/samplesheet.csv --outdir results
```

Full step-by-step guide, including Seqera Cloud launch: [docs/NEXTFLOW_TUTORIAL.md](docs/NEXTFLOW_TUTORIAL.md).

## Requirements

- Nextflow >= 24.10
- A container engine: Apptainer (HPC/Slurm) or Docker
- Java 17+ (Nextflow runtime)
