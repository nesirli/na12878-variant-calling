/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: Variant filtering and merging
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Mirrors the Snakemake rules:
      SelectVariants (SNP/INDEL) -> VariantFiltration -> keep PASS (bcftools view)
      -> bcftools concat -> bcftools stats.
----------------------------------------------------------------------------------------
*/

include { GATK4_SELECTVARIANTS   as GATK4_SELECTVARIANTS_SNP    } from '../../modules/nf-core/gatk4/selectvariants'
include { GATK4_SELECTVARIANTS   as GATK4_SELECTVARIANTS_INDEL  } from '../../modules/nf-core/gatk4/selectvariants'
include { GATK4_VARIANTFILTRATION as GATK4_VARIANTFILTRATION_SNP   } from '../../modules/nf-core/gatk4/variantfiltration'
include { GATK4_VARIANTFILTRATION as GATK4_VARIANTFILTRATION_INDEL } from '../../modules/nf-core/gatk4/variantfiltration'
include { BCFTOOLS_VIEW    as BCFTOOLS_VIEW_SNP    } from '../../modules/nf-core/bcftools/view'
include { BCFTOOLS_VIEW    as BCFTOOLS_VIEW_INDEL  } from '../../modules/nf-core/bcftools/view'
include { BCFTOOLS_CONCAT  } from '../../modules/nf-core/bcftools/concat'
include { BCFTOOLS_STATS   } from '../../modules/nf-core/bcftools/stats'

workflow FILTER_VARIANTS {

    take:
    ch_vcf                      // channel: [ meta, raw.vcf.gz ]
    ch_vcf_index                // channel: [ meta, raw.vcf.gz.tbi ]
    ch_fasta                    // channel: [ meta, fasta ]
    ch_fai                      // channel: [ meta, fai ]
    ch_dict                     // channel: [ meta, dict ]

    main:
    //
    // MODULE: Select SNPs / INDELs
    //
    def ch_vcf_idx = ch_vcf.join(ch_vcf_index).map { meta, vcf, tbi -> [ meta, vcf, tbi, [] ] }

    GATK4_SELECTVARIANTS_SNP(ch_vcf_idx)
    GATK4_SELECTVARIANTS_INDEL(ch_vcf_idx)

    //
    // MODULE: VariantFiltration
    //
    GATK4_VARIANTFILTRATION_SNP(
        GATK4_SELECTVARIANTS_SNP.out.vcf.join(GATK4_SELECTVARIANTS_SNP.out.tbi),
        ch_fasta, ch_fai, ch_dict, [ [], [] ]
    )
    GATK4_VARIANTFILTRATION_INDEL(
        GATK4_SELECTVARIANTS_INDEL.out.vcf.join(GATK4_SELECTVARIANTS_INDEL.out.tbi),
        ch_fasta, ch_fai, ch_dict, [ [], [] ]
    )

    //
    // MODULE: Keep only PASS variants
    //
    BCFTOOLS_VIEW_SNP(GATK4_VARIANTFILTRATION_SNP.out.vcf.join(GATK4_VARIANTFILTRATION_SNP.out.tbi), [], [], [])
    BCFTOOLS_VIEW_INDEL(GATK4_VARIANTFILTRATION_INDEL.out.vcf.join(GATK4_VARIANTFILTRATION_INDEL.out.tbi), [], [], [])

    //
    // MODULE: Merge SNPs and indels
    //
    def ch_concat_in = BCFTOOLS_VIEW_SNP.out.vcf
        .combine(BCFTOOLS_VIEW_SNP.out.index)
        .combine(BCFTOOLS_VIEW_INDEL.out.vcf)
        .combine(BCFTOOLS_VIEW_INDEL.out.index)
        .map { m1, v1, _m2, t1, _m3, v2, _m4, t2 ->
            [ [ id: m1.id ], [ v1, v2 ], [ t1, t2 ] ]
        }
    BCFTOOLS_CONCAT(ch_concat_in)

    //
    // MODULE: Variant statistics (replaces the SNP/indel/TsTv text files)
    //
    BCFTOOLS_STATS(
        BCFTOOLS_CONCAT.out.vcf.join(BCFTOOLS_CONCAT.out.index).map { meta, vcf, tbi -> [ meta, vcf, tbi ] },
        [ [], [] ], [ [], [] ], [ [], [] ], [ [], [] ], [ [], [] ]
    )

    emit:
    vcf     = BCFTOOLS_CONCAT.out.vcf                // channel: [ meta, final.vcf.gz ]
    vcf_idx = BCFTOOLS_CONCAT.out.index              // channel: [ meta, final.vcf.gz.tbi ]
    stats   = BCFTOOLS_STATS.out.stats               // channel: [ meta, bcftools_stats.txt ]
}
