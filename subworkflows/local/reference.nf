/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: Reference preparation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Obtains the GRCh38 primary assembly, the BQSR known sites, and prepares the
    indexes the alignment and GATK stages need:
      - samtools faidx  -> .fai (+ .sizes)
      - bwa index       -> .amb/.ann/.bwt/.pac/.sa
----------------------------------------------------------------------------------------
*/

include { DOWNLOAD_REFERENCE } from '../../modules/local/download_reference'
include { SAMTOOLS_FAIDX     } from '../../modules/nf-core/samtools/faidx'
include { BWA_INDEX          } from '../../modules/nf-core/bwa/index'

workflow PREPARE_REFERENCE {

    take:
    ch_reference                // channel: [ meta ]

    main:
    DOWNLOAD_REFERENCE(ch_reference)

    def ch_fasta = DOWNLOAD_REFERENCE.out.fasta

    //
    // MODULE: Create FASTA index
    //
    SAMTOOLS_FAIDX(ch_fasta.map { meta, fasta -> [ meta, fasta, [] ] }, true)

    //
    // MODULE: Build BWA index
    //
    BWA_INDEX(ch_fasta)

    emit:
    fasta      = ch_fasta                            // channel: [ meta, fasta ]
    fai        = SAMTOOLS_FAIDX.out.fai              // channel: [ meta, fai ]
    sizes      = SAMTOOLS_FAIDX.out.sizes            // channel: [ meta, sizes ]
    bwa_index  = BWA_INDEX.out.index                 // channel: [ meta, index_dir ]
    dbsnp      = DOWNLOAD_REFERENCE.out.dbsnp        // channel: [ meta, dbsnp.vcf, dbsnp.vcf.idx ]
    indels     = DOWNLOAD_REFERENCE.out.indels       // channel: [ meta, known_indels.vcf.gz, known_indels.vcf.gz.tbi ]
    integrity  = DOWNLOAD_REFERENCE.out.integrity    // path: reference_integrity.txt
}
