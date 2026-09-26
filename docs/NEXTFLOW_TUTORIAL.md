# Nextflow + nf-core port: step-by-step tutorial

This document teaches how the `na12878-variant-calling` analysis was ported from
Snakemake to Nextflow using nf-core modules, and how to run it on a Slurm cluster
launched from Seqera Cloud.

> Work in progress: each section is added as the corresponding stage is implemented.

## Contents

1. [Why nf-core modules](#1-why-nf-core-modules)
2. [Pipeline skeleton](#2-pipeline-skeleton)
3. [Read QC](#3-read-qc)
4. [Alignment](#4-alignment)
5. [GATK variant calling](#5-gatk-variant-calling)
6. [Filtering and merging](#6-filtering-and-merging)
7. [Annotation](#7-annotation)
8. [Validation against GIAB](#8-validation-against-giab)
9. [Execution profiles and resources](#9-execution-profiles-and-resources)
10. [Running on Slurm](#10-running-on-slurm)
11. [Launching from Seqera Cloud](#11-launching-from-seqera-cloud)

## 1. Why nf-core modules

An nf-core module is a self-contained Nextflow process for one tool (or one tool
sub-command). Each module ships:

- `main.nf` — the process definition with typed `input:`/`output:` and a `meta` map.
- `environment.yml` and a container reference (BioContainers) for reproducibility.
- a `versions.yml` emitter so tool versions are tracked per run.
- `<tool>.config` describing resource hints.

Reusing them means we inherit tested command lines, container images, and version
tracking instead of reimplementing wrapper scripts. We install them with
`nf-core modules install <tool/subcommand>` and call them from a normal Nextflow
workflow, passing extra CLI flags through `ext.args` in a config file.

## 3. Read QC

The first real stage is the simplest: FastQC per sample, then MultiQC. It teaches the
three moving parts of every stage: install a module, call it from a subworkflow, and
feed its reports to MultiQC.

Install the module (from the pipeline root):

```bash
nf-core modules install fastqc
```

`FASTQC` takes `[ meta, [ fastq_1, fastq_2 ] ]` — the sample metadata map plus the list
of FASTQ paths, which is exactly what `PIPELINE_INITIALISATION` produces from the
samplesheet. It emits `html`, `zip`, and (via the version topic) `versions_fastqc`.

`subworkflows/local/qc.nf` wraps it:

```groovy
include { FASTQC } from '../../modules/nf-core/fastqc/main'

workflow QC {
    take:
    ch_samplesheet

    main:
    def ch_multiqc_files = channel.empty()
    FASTQC(ch_samplesheet)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map { _meta, zip -> zip })

    emit:
    multiqc_files = ch_multiqc_files
}
```

Two nf-core conventions worth noting:

- **Versions are collected through a channel topic**, not by wiring `versions` outputs
  by hand. Modules publish `emit: ..., topic: versions`, and the top-level workflow
  reads `channel.topic("versions")`. That is why the subworkflow above emits only the
  MultiQC inputs.
- **`meta` propagates everywhere.** Every module output is `[ meta, files... ]` so
  sample identity survives the whole DAG without global variables.

`MULTIQC` stays in the top-level workflow so it can aggregate reports from *all* stages
(FastQC now, alignment/variant metrics later).

Test it without running any tools using a stub run:

```bash
nextflow run . -stub-run --input assets/samplesheet.csv --outdir results
```

Commit: `feat(qc): add FastQC subworkflow`.

<!-- sections below are filled in as stages land -->

