/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: Functional annotation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SnpEff annotation plus an awk-based extraction of high-impact and missense
    variants (the bioconda SnpEff package does not ship SnpSift).
----------------------------------------------------------------------------------------
*/

include { SNPEFF_DOWNLOAD } from '../../modules/nf-core/snpeff/download'
include { SNPEFF_SNPEFF   } from '../../modules/nf-core/snpeff/snpeff'
include { ANN_FILTER      } from '../../modules/local/ann_filter'

workflow ANNOTATION {

    take:
    ch_vcf                      // channel: [ meta, final.vcf.gz ]

    main:
    def ch_snpeff_db = channel.of([ [ id: 'snpeff' ], params.snpeff_db ])

    //
    // MODULE: Download the SnpEff database
    //
    SNPEFF_DOWNLOAD(ch_snpeff_db)

    //
    // MODULE: Annotate variants
    //
    SNPEFF_SNPEFF(
        ch_vcf,
        params.snpeff_db,
        SNPEFF_DOWNLOAD.out.cache
    )

    //
    // MODULE: Extract high-impact and missense variants
    //
    ANN_FILTER(
        SNPEFF_SNPEFF.out.vcf,
        file("${projectDir}/scripts/filter_ann.awk", checkIfExists: true)
    )

    emit:
    vcf         = SNPEFF_SNPEFF.out.vcf            // channel: [ meta, annotated.vcf ]
    high_impact = ANN_FILTER.out.high_impact       // channel: [ meta, high_impact.vcf ]
    missense    = ANN_FILTER.out.missense          // channel: [ meta, missense.vcf ]
    summary     = ANN_FILTER.out.summary           // channel: [ meta, annotation_summary.txt ]
    report      = SNPEFF_SNPEFF.out.report         // channel: multiqc files
    summary_html = SNPEFF_SNPEFF.out.summary_html  // channel: multiqc files
    genes_txt   = SNPEFF_SNPEFF.out.genes_txt      // channel: multiqc files
}
