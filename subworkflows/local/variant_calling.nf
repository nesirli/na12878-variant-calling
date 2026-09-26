/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: GATK variant calling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Mirrors the Snakemake analysis rules:
      CreateSequenceDictionary -> MarkDuplicates -> BaseRecalibrator -> ApplyBQSR
      -> HaplotypeCaller.
    The Broad known-sites VCFs use UCSC contig names (chr20); the Ensembl reference
    uses (20). They are renamed with bcftools/annotate before BQSR.
----------------------------------------------------------------------------------------
*/

include { GATK4_CREATESEQUENCEDICTIONARY } from '../../modules/nf-core/gatk4/createsequencedictionary'
include { GATK4_MARKDUPLICATES            } from '../../modules/nf-core/gatk4/markduplicates'
include { GATK4_BASERECALIBRATOR          } from '../../modules/nf-core/gatk4/baserecalibrator'
include { GATK4_APPLYBQSR                 } from '../../modules/nf-core/gatk4/applybqsr'
include { GATK4_HAPLOTYPECALLER           } from '../../modules/nf-core/gatk4/haplotypecaller'
include { BCFTOOLS_ANNOTATE as BCFTOOLS_ANNOTATE_DBSNP  } from '../../modules/nf-core/bcftools/annotate'
include { BCFTOOLS_ANNOTATE as BCFTOOLS_ANNOTATE_INDELS } from '../../modules/nf-core/bcftools/annotate'
include { BCFTOOLS_INDEX    as BCFTOOLS_INDEX_DBSNP     } from '../../modules/nf-core/bcftools/index'
include { BCFTOOLS_INDEX    as BCFTOOLS_INDEX_INDELS    } from '../../modules/nf-core/bcftools/index'

workflow VARIANT_CALLING {

    take:
    ch_bam                      // channel: [ meta, bam ]
    ch_fasta                    // channel: [ meta, fasta ]
    ch_fai                      // channel: [ meta, fai ]
    ch_dbsnp                    // channel: [ meta, dbsnp.vcf, dbsnp.vcf.idx ]
    ch_indels                   // channel: [ meta, known_indels.vcf.gz, known_indels.vcf.gz.tbi ]

    main:
    //
    // MODULE: Sequence dictionary
    //
    GATK4_CREATESEQUENCEDICTIONARY(ch_fasta)
    def ch_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict

    //
    // MODULE: Rename known sites to Ensembl contig names
    //
    def ch_chr_map = file("${projectDir}/assets/chr_rename.txt", checkIfExists: true)

    def ch_dbsnp_in = ch_dbsnp.map { _meta, vcf, idx ->
        [ [ id: 'Homo_sapiens_assembly38.dbsnp138.ensembl' ], vcf, idx, [], [], [], [], ch_chr_map ]
    }
    def ch_indels_in = ch_indels.map { _meta, vcf, tbi ->
        [ [ id: 'Homo_sapiens_assembly38.known_indels.ensembl' ], vcf, tbi, [], [], [], [], ch_chr_map ]
    }

    BCFTOOLS_ANNOTATE_DBSNP(ch_dbsnp_in)
    def ch_dbsnp_renamed = BCFTOOLS_ANNOTATE_DBSNP.out.vcf
    BCFTOOLS_INDEX_DBSNP(ch_dbsnp_renamed)
    def ch_dbsnp_renamed_idx = BCFTOOLS_INDEX_DBSNP.out.index

    BCFTOOLS_ANNOTATE_INDELS(ch_indels_in)
    def ch_indels_renamed = BCFTOOLS_ANNOTATE_INDELS.out.vcf
    BCFTOOLS_INDEX_INDELS(ch_indels_renamed)
    def ch_indels_renamed_idx = BCFTOOLS_INDEX_INDELS.out.index

    // Both known-sites files presented to BaseRecalibrator as lists
    def ch_known_sites = ch_dbsnp_renamed.combine(ch_indels_renamed)
        .map { _m1, v1, _m2, v2 -> [ [ id: 'known_sites' ], [ v1, v2 ] ] }
    def ch_known_sites_tbi = ch_dbsnp_renamed_idx.combine(ch_indels_renamed_idx)
        .map { _m1, t1, _m2, t2 -> [ [ id: 'known_sites' ], [ t1, t2 ] ] }

    //
    // MODULE: Mark duplicates
    //
    GATK4_MARKDUPLICATES(ch_bam, ch_fasta.map { _m, fa -> fa }, ch_fai.map { _m, fai -> fai })
    def ch_dedup_bam = GATK4_MARKDUPLICATES.out.bam
    def ch_dedup_bai = GATK4_MARKDUPLICATES.out.bai

    //
    // MODULE: BaseRecalibrator
    //
    GATK4_BASERECALIBRATOR(
        ch_dedup_bam.join(ch_dedup_bai).map { meta, bam, bai -> [ meta, bam, bai, [] ] },
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_known_sites,
        ch_known_sites_tbi
    )
    def ch_recal_table = GATK4_BASERECALIBRATOR.out.table

    //
    // MODULE: ApplyBQSR
    //
    GATK4_APPLYBQSR(
        ch_dedup_bam.join(ch_dedup_bai).join(ch_recal_table).map { meta, bam, bai, table -> [ meta, bam, bai, table, [] ] },
        ch_fasta.join(ch_fai).join(ch_dict).map { meta, fa, fai, dict -> [ meta, fa, fai, dict ] },
        'bam'
    )
    def ch_recal_bam = GATK4_APPLYBQSR.out.bam
    def ch_recal_bai = GATK4_APPLYBQSR.out.bai

    //
    // MODULE: HaplotypeCaller
    //
    GATK4_HAPLOTYPECALLER(
        ch_recal_bam.join(ch_recal_bai).map { meta, bam, bai -> [ meta, bam, bai, [], [] ] },
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_dbsnp_renamed,
        ch_dbsnp_renamed_idx
    )

    emit:
    dict        = ch_dict                                   // channel: [ meta, dict ]
    dedup_bam   = ch_dedup_bam                              // channel: [ meta, bam ]
    dedup_bai   = ch_dedup_bai                              // channel: [ meta, bai ]
    recal_table = ch_recal_table                            // channel: [ meta, table ]
    recal_bam   = ch_recal_bam                              // channel: [ meta, bam ]
    recal_bai   = ch_recal_bai                              // channel: [ meta, bai ]
    vcf         = GATK4_HAPLOTYPECALLER.out.vcf             // channel: [ meta, vcf.gz ]
    vcf_index   = GATK4_HAPLOTYPECALLER.out.tbi             // channel: [ meta, tbi ]
}
