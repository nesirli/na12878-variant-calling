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

## 4. Alignment

Alignment needs a reference and its indexes first. We download GRCh38 in a **local
module** (`modules/local/download_reference`) because no nf-core module covers the
Ensembl FASTA + BSD-`sum` verification + Broad known-sites combination. Writing a local
module is the same contract as an nf-core module: typed `input:`, `output:` with
`emit:`, a container, and a `versions.yml`.

`subworkflows/local/reference.nf` then uses two nf-core modules:

```groovy
SAMTOOLS_FAIDX(ch_fasta.map { meta, fasta -> [ meta, fasta, [] ] }, true)  // -> .fai (+ .sizes)
BWA_INDEX(ch_fasta)                                                        // -> index dir
```

`subworkflows/local/align.nf` mirrors the Snakemake rules:

```groovy
BWA_MEM(ch_reads, ch_bwa_index, ch_fasta, false)   // false = do not sort inside BWA
SAMTOOLS_SORT(BWA_MEM.out.bam, [ [:], [], [] ], '')
SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)
SAMTOOLS_FLAGSTAT(SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.index))
```

Four things are worth internalising here:

- **`meta` is the join key.** `ch_reads` carries the sample `meta`; `ch_bwa_index` and
  `ch_fasta` carry a reference `meta`. Nextflow pairs each sample with the single
  reference entry automatically. `join` pairs the BAM with its index for `flagstat`.
- **`ext.args` carries tool flags.** The BWA read group is injected in
  `conf/modules.config` rather than hard-coded in the module:
  `ext.args = { "-R '@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:lib1'" }`.
- **`ext.prefix` controls output names.** `SAMTOOLS_SORT` is told to write
  `${meta.id}.sorted.bam` so it never collides with the unsorted `BWA_MEM` output.
- **Resource labels, not hard numbers.** Modules declare `label 'process_low'`,
  `process_high`, etc.; `conf/base.config` maps those to cpus/memory/time. We lowered
  the template's defaults (72/200 GB) to values that fit this 47 GB cluster, and added
  `resourceLimits` for tests.

Stub run (validates the whole DAG without running tools):

```bash
nextflow run . -stub-run -c conf/local_test.config --input assets/samplesheet.csv --outdir results
```

Commit: `feat(align): add reference prep and BWA-MEM alignment subworkflows`.

## 5. GATK variant calling

`subworkflows/local/variant_calling.nf` chains the GATK4 modules the same way the
Snakemake `04_analysis.smk` rules do:

```groovy
GATK4_CREATESEQUENCEDICTIONARY(ch_fasta)                 // -> .dict
GATK4_MARKDUPLICATES(ch_bam, fasta, fai)                 // -> dedup.bam/.bai/.metrics
GATK4_BASERECALIBRATOR(dedup + index, fasta, fai, dict, known_sites, known_tbi)
GATK4_APPLYBQSR(dedup + index + table, fasta + fai + dict, 'bam')
GATK4_HAPLOTYPECALLER(recal + index, fasta, fai, dict, dbsnp, dbsnp_tbi)
```

Lessons this stage teaches:

- **Aliased includes for repeated tools.** Nextflow refuses to use the same process
  twice in one scope (`Process 'BCFTOOLS_ANNOTATE' has been already used`). Rename the
  known-sites with two aliases:
  `include { BCFTOOLS_ANNOTATE as BCFTOOLS_ANNOTATE_DBSNP }` and `..._INDELS`. Config
  selectors then use a regex: `withName: 'BCFTOOLS_ANNOTATE.*'`.
- **Contig harmonisation.** The Broad known-sites VCFs use `chr20`, the Ensembl reference
  uses `20`. `assets/chr_rename.txt` is passed to `bcftools/annotate --rename-chrs`; the
  module exposes this as the `rename_chrs` input slot.
- **Avoiding a module quirk.** `bcftools/annotate` runs `bcftools index <input>` when its
  `index` input is empty. We pass the existing `.idx`/`.tbi` so that pre-command is
  skipped, then index the *output* with a separate `BCFTOOLS_INDEX_*` call.
- **`ext.args` is where pipeline policy lives.** All the GATK flags are in
  `conf/modules.config`:
  `GATK4_MARKDUPLICATES` gets `--REMOVE_DUPLICATES false --CREATE_INDEX true`;
  `GATK4_HAPLOTYPECALLER` gets `-L ${params.variant_interval}` (chr20 by default).
- **Multi-file inputs are lists.** `BaseRecalibrator` takes a list of known sites, so the
  two renamed VCFs are combined with `combine()` into `[ [id:'known_sites'], [v1,v2] ]`.

Commit: `feat(variant-calling): add GATK4 MarkDuplicates/BQSR/HaplotypeCaller`.

## 6. Filtering and merging

`subworkflows/local/filter_variants.nf`:

```groovy
GATK4_SELECTVARIANTS_SNP(...)                              // --select-type-to-include SNP
GATK4_VARIANTFILTRATION_SNP(...)                           // QD/MQ/FS/SOR filters
BCFTOOLS_VIEW_SNP(filtered, [], [], [])                    // -f PASS
GATK4_SELECTVARIANTS_INDEL(...) / GATK4_VARIANTFILTRATION_INDEL(...)
BCFTOOLS_CONCAT(...)                                        // -a -Oz, merge SNP + INDEL
BCFTOOLS_STATS(...)
```

Notes:

- **Empty tuple inputs must have the right arity.** `VariantFiltration` declares
  `tuple val(meta), path(gzi)` for the optional `.gzi`; passing `[]` fails with
  `Input tuple does not match tuple declaration`. Pass `[ [], [] ]` instead. The same
  applies to `bcftools/stats`' optional `[meta, file]` inputs.
- **`--write-index` keeps indexes tied to their VCFs.** `bcftools view` and `concat` are
  configured with `--write-index=tbi`, so each emitted VCF always arrives with its `.tbi`
  in the same tuple — no dangling index between processes.
- **One stats file replaces three text files.** `bcftools stats` reports SNP/indel counts
  and the Ts/Tv ratio that the Snakemake rules wrote separately.

Commit: `feat(filter): add GATK SelectVariants/VariantFiltration + bcftools merge`.

## 7. Annotation

`subworkflows/local/annotation.nf`:

```groovy
SNPEFF_DOWNLOAD(ch_snpeff_db)                 // [ [id:'snpeff'], params.snpeff_db ]
SNPEFF_SNPEFF(ch_vcf, params.snpeff_db, SNPEFF_DOWNLOAD.out.cache)
ANN_FILTER(SNPEFF_SNPEFF.out.vcf, script)     // local module
```

- **Database as a first-class input.** The SnpEff database name (`GRCh38.105`) is a
  param. `snpeff/download` emits a cache directory which is passed as `meta2` to
  `snpeff/snpeff`; Nextflow stages it into the task working directory and the module
  points `-dataDir` at it.
- **Why a local module here.** The bioconda SnpEff package does not ship SnpSift, so
  selecting high-impact/missense variants is done with `scripts/filter_ann.awk`. Wrapping
  it as `modules/local/ann_filter` keeps the same container/version/reporting contract as
  an nf-core module.
- **MultiQC via a channel topic.** SnpEff publishes its `.csv`/`.html`/`.genes.txt` to the
  `multiqc_files` topic. The top-level workflow subscribes with
  `channel.topic("multiqc_files")` and mixes the file paths into the MultiQC inputs — no
  manual wiring between the annotation subworkflow and MultiQC.

Commit: `feat(annotate): add SnpEff annotation and impact filtering`.

<!-- sections below are filled in as stages land -->





