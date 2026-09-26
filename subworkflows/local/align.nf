/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: Read alignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BWA-MEM against the GRCh38 reference, then coordinate-sort, index, and
    alignment summary statistics with samtools.
----------------------------------------------------------------------------------------
*/

include { BWA_MEM           } from '../../modules/nf-core/bwa/mem'
include { SAMTOOLS_SORT     } from '../../modules/nf-core/samtools/sort'
include { SAMTOOLS_INDEX    } from '../../modules/nf-core/samtools/index'
include { SAMTOOLS_FLAGSTAT } from '../../modules/nf-core/samtools/flagstat'

workflow ALIGN {

    take:
    ch_reads                    // channel: [ meta, [ fastq_1, fastq_2 ] ]
    ch_fasta                    // channel: [ meta_ref, fasta ]
    ch_bwa_index                // channel: [ meta_ref, index_dir ]

    main:
    //
    // MODULE: BWA-MEM (unsorted BAM via `samtools view`)
    //
    BWA_MEM(ch_reads, ch_bwa_index, ch_fasta, false)
    def ch_bam = BWA_MEM.out.bam

    //
    // MODULE: Coordinate sort
    //
    SAMTOOLS_SORT(ch_bam, [ [:], [], [] ], '')

    //
    // MODULE: Index the sorted BAM
    //
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)

    //
    // MODULE: flagstat summary
    //
    SAMTOOLS_FLAGSTAT(SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.index))

    emit:
    bam    = SAMTOOLS_SORT.out.bam          // channel: [ meta, bam ]
    bai    = SAMTOOLS_INDEX.out.index       // channel: [ meta, bai ]
    stats  = SAMTOOLS_FLAGSTAT.out.flagstat // channel: [ meta, flagstat ]
}
