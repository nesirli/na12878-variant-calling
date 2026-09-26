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

<!-- sections below are filled in as stages land -->
