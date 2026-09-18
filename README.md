# na12878-variant-calling

An end-to-end reproducible workflow for calling, filtering, and annotating clinical genomic variants.

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/) >= 9
- Either a conda-based solver (environments are pinned in `envs/`) or [Apptainer](https://apptainer.org/) for containerized execution

## Configuration

Edit `config/config.yaml` to set input/output directories and per-rule thread and memory parameters.

## Usage

Each rule declares both a `conda:` environment and a BioContainers `container:` image, so pick one backend per run (do not pass both to `--sdm`).

Run the workflow with conda-managed environments on all available cores:

```bash
snakemake --cores 12 --sdm conda \
  --resources mem_mb=45000 \
  --scheduler ilp \
  --rerun-incomplete --keep-going --retries 3 \
  --latency-wait 60
```

Run the workflow with Apptainer (no conda environments required):

```bash
snakemake --cores 12 --sdm apptainer \
  --resources mem_mb=45000 \
  --scheduler ilp \
  --rerun-incomplete --keep-going --retries 3 \
  --latency-wait 60
```

Dry-run to preview the DAG:

```bash
snakemake --dryrun --cores 12 --sdm conda --resources mem_mb=45000 --scheduler ilp
```

Hint: on HPC systems where `$HOME` is small or read-only, point Apptainer at a large scratch location, e.g. `--apptainer-prefix /scratch/$USER/snakemake-apptainer`.

Note: the BioContainers images provide `awk`, `gzip`, `gunzip`, `sum`, `tee`, and `mkdir` via BusyBox, while the conda backend pins the GNU equivalents (`coreutils=9.4`, `gzip=1.14`, `gawk=5.4.1`). The operations used here are functionally equivalent — the Ensembl `sum` checksum was verified to produce identical output on both.

## Outputs


- `data/raw/` — downloaded reads and reference
- `results/qc/` — FastQC and MultiQC reports
- `results/align/` — sorted BAM files, indexes, and alignment stats
- `logs/` — per-rule logs
