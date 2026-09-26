/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: Reference preparation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Obtains the GRCh38 primary assembly and prepares the indexes the alignment
    and GATK stages need:
      - samtools faidx  -> .fai
      - bwa index       -> .amb/.ann/.bwt/.pac/.sa
      - (GATK .dict is added in the variant-calling subworkflow)
----------------------------------------------------------------------------------------
*/

include { DOWNLOAD_REFERENCE } from '../../modules/local/download_reference'
include { SAMTOOLS_FAIDX     } from '../../modules/nf-core/samtools/faidx'
include { BWA_INDEX          } from '../../modules/nf-core/bwa/index'

workflow PREPARE_REFERENCE {

    take:
    ch_reference                // channel: [ meta ]

    main:
    def ch_fasta       = channel.empty()
    def ch_known_sites = channel.empty()

    if (params.download_reference) {
        DOWNLOAD_REFERENCE(ch_reference)
        ch_fasta       = DOWNLOAD_REFERENCE.out.fasta
        ch_known_sites = DOWNLOAD_REFERENCE.out.known_sites
    } else {
        ch_fasta = ch_reference.map { meta ->
            [ meta, file(params.fasta, checkIfExists: true) ]
        }
    }

    //
    // MODULE: Create FASTA index
    //
    SAMTOOLS_FAIDX(ch_fasta.map { meta, fasta -> [ meta, fasta, [] ] }, true)

    //
    // MODULE: Build BWA index
    //
    BWA_INDEX(ch_fasta)

    emit:
    fasta       = ch_fasta                                  // channel: [ meta, fasta ]
    fai         = SAMTOOLS_FAIDX.out.fai                    // channel: [ meta, fai ]
    sizes       = SAMTOOLS_FAIDX.out.sizes                  // channel: [ meta, sizes ]
    bwa_index   = BWA_INDEX.out.index                       // channel: [ meta, index_dir ]
    known_sites = ch_known_sites                            // channel: [ meta, [ vcf... ] ]
}
