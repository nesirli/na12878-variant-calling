/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: Validation against the GIAB truth set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Downloads the NA12878/GH001 (HG001) GIAB benchmark, renames its UCSC contigs to
    Ensembl, intersects with the calls (bcftools isec), and reports sensitivity and
    precision.
----------------------------------------------------------------------------------------
*/

include { WGET              as WGET_TRUTH_VCF } from '../../modules/nf-core/wget'
include { WGET              as WGET_TRUTH_TBI } from '../../modules/nf-core/wget'
include { BCFTOOLS_ANNOTATE as BCFTOOLS_ANNOTATE_TRUTH } from '../../modules/nf-core/bcftools/annotate'
include { BCFTOOLS_INDEX                                   } from '../../modules/nf-core/bcftools/index'
include { BCFTOOLS_ISEC                                    } from '../../modules/nf-core/bcftools/isec'
include { VALIDATE_METRICS  } from '../../modules/local/validate_metrics'

workflow VALIDATION {

    take:
    ch_vcf                      // channel: [ meta, final.vcf.gz ]
    ch_vcf_index                // channel: [ meta, final.vcf.gz.tbi ]

    main:
    //
    // MODULE: Download the GIAB truth set (VCF + index)
    //
    def ch_truth_vcf = channel.of([
        [ id: params.truth_filename.replaceAll(/\.vcf\.gz$/, '') ],
        "${params.truth_url}/${params.truth_filename}",
        'vcf.gz'
    ])
    def ch_truth_tbi = channel.of([
        [ id: "${params.truth_filename}" ],
        "${params.truth_url}/${params.truth_filename}.tbi",
        'tbi'
    ])
    WGET_TRUTH_VCF(ch_truth_vcf)
    WGET_TRUTH_TBI(ch_truth_tbi)

    //
    // MODULE: Rename truth-set contigs to Ensembl style
    //
    def ch_chr_map = file("${projectDir}/assets/chr_rename.txt", checkIfExists: true)
    def ch_truth_in = WGET_TRUTH_VCF.out.outfile
        .combine(WGET_TRUTH_TBI.out.outfile)
        .map { _m1, vcf, _m2, tbi ->
            [ [ id: 'truth.ensembl' ], vcf, tbi, [], [], [], [], ch_chr_map ]
        }
    BCFTOOLS_ANNOTATE_TRUTH(ch_truth_in)
    BCFTOOLS_INDEX(BCFTOOLS_ANNOTATE_TRUTH.out.vcf)
    def ch_truth_renamed     = BCFTOOLS_ANNOTATE_TRUTH.out.vcf
    def ch_truth_renamed_idx = BCFTOOLS_INDEX.out.index

    //
    // MODULE: Intersect calls and truth (bcftools isec)
    //
    def ch_isec_in = ch_vcf.join(ch_vcf_index)
        .combine(ch_truth_renamed.join(ch_truth_renamed_idx))
        .map { meta, calls, calls_idx, _tm, truth, truth_idx ->
            [ meta, [ calls, truth ], [ calls_idx, truth_idx ], [], [], [] ]
        }
    BCFTOOLS_ISEC(ch_isec_in)

    //
    // MODULE: Concordance metrics
    //
    VALIDATE_METRICS(BCFTOOLS_ISEC.out.results)

    emit:
    report       = VALIDATE_METRICS.out.report       // channel: [ meta, validation.txt ]
    isec_results = BCFTOOLS_ISEC.out.results         // channel: [ meta, isec_dir ]
}
