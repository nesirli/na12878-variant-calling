/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { DOWNLOAD_READS         } from '../modules/local/download_reads'
include { QC                     } from '../subworkflows/local/qc'
include { PREPARE_REFERENCE      } from '../subworkflows/local/reference'
include { ALIGN                  } from '../subworkflows/local/align'
include { VARIANT_CALLING        } from '../subworkflows/local/variant_calling'
include { FILTER_VARIANTS        } from '../subworkflows/local/filter_variants'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_na12878_variant_calling_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NA12878_VARIANT_CALLING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // Reads: download from ENA, or use FASTQ paths listed in the samplesheet
    //
    def ch_reads = channel.empty()
    if (params.download_reads) {
        DOWNLOAD_READS(ch_samplesheet.map { meta, _reads -> meta })
        ch_reads = DOWNLOAD_READS.out.reads
    } else {
        ch_reads = ch_samplesheet
    }

    //
    // SUBWORKFLOW: Prepare the reference indexes (faidx + bwa index)
    //
    def ch_reference = channel.of([ id: params.reference_name ?: 'GRCh38' ])
    PREPARE_REFERENCE(ch_reference)

    //
    // SUBWORKFLOW: Read QC (FastQC)
    //
    QC(ch_reads)
    ch_multiqc_files = ch_multiqc_files.mix(QC.out.multiqc_files)

    //
    // SUBWORKFLOW: Read alignment (BWA-MEM + samtools)
    //
    ALIGN(ch_reads, PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.bwa_index)

    //
    // SUBWORKFLOW: GATK variant calling (MarkDuplicates, BQSR, HaplotypeCaller)
    //
    VARIANT_CALLING(
        ALIGN.out.bam,
        PREPARE_REFERENCE.out.fasta,
        PREPARE_REFERENCE.out.fai,
        PREPARE_REFERENCE.out.dbsnp,
        PREPARE_REFERENCE.out.indels
    )

    //
    // SUBWORKFLOW: Filter and merge variants (GATK + bcftools)
    //
    FILTER_VARIANTS(
        VARIANT_CALLING.out.vcf,
        VARIANT_CALLING.out.vcf_index,
        PREPARE_REFERENCE.out.fasta,
        PREPARE_REFERENCE.out.fai,
        VARIANT_CALLING.out.dict
    )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'na12878-variant-calling_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'na12878-variant-calling'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
