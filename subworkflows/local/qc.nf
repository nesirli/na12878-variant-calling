/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: Read quality control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Runs FastQC over the raw reads of every sample and collects the per-sample
    reports so the top-level workflow can feed them to MultiQC.
----------------------------------------------------------------------------------------
*/

include { FASTQC } from '../../modules/nf-core/fastqc/main'

workflow QC {

    take:
    ch_samplesheet              // channel: [ meta, [ fastq_1, fastq_2 ] ]

    main:
    def ch_multiqc_files = channel.empty()

    //
    // MODULE: Run FastQC
    //
    FASTQC(ch_samplesheet)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map { _meta, zip -> zip })

    emit:
    multiqc_files = ch_multiqc_files  // channel: [ path(zip) ]
}
