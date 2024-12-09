include { ACI                                               } from '../../modules/local/aci'  
include { BCFTOOLS                                          } from '../../modules/local/bcftools'
include { IGV_REPORTS                                       } from '../../modules/local/igvreports' 
include { IVAR_VARIANTS                                     } from '../../modules/local/ivar'    
include { FASTQC                                            } from '../../modules/local/fastqc' 
include { KRAKEN2                                           } from '../../modules/local/kraken2' 
include { SAMTOOLS_QC as INITIAL_QC                         } from '../../modules/local/samtools'  
include { SAMTOOLS_QC as FINAL_QC                           } from '../../modules/local/samtools'   
include { SAMTOOLS_AMPLICONSTATS as AMPLICONSTATS           } from '../../modules/local/samtools'
include { SAMTOOLS_PLOT_AMPLICONSTATS as PLOT_AMPLICONSTATS } from '../../modules/local/samtools' 


workflow QC {
  take:
    ch_raw_reads // channel: [meta, reads]
    ch_clean_reads // channel: [meta, reads]
    ch_kraken2_db // channel: path
    ch_initial_bam // channel: [meta, sam]
    ch_trim_bam // channel: [meta, bam, bai]
    ch_reference_genome // channel: fasta
    ch_gff_file // channel: gff file
    ch_amplicon_bed  // channel: bedfile
    ch_primer_bed  // channel: bedfile

  main:
    ch_for_multiqc = Channel.empty()
    ch_summary     = Channel.empty()
    ch_versions    = Channel.empty()
    
    if (params.fastqc ) {
      FASTQC(ch_raw_reads)
      ch_summary     = ch_summary.mix(FASTQC.out.fastq_name.collectFile(name: "fastq_names.csv", keepHeader: true ))
      ch_versions    = ch_versions.mix(FASTQC.out.versions.first())
      ch_for_multiqc = ch_for_multiqc.mix(FASTQC.out.fastqc_files)
    }

    if ( params.kraken2_db ) {
      KRAKEN2(ch_clean_reads.combine(ch_kraken2_db))
      ch_for_multiqc = ch_for_multiqc.mix(KRAKEN2.out.kraken2_files)
      ch_summary     = ch_summary.mix(KRAKEN2.out.kraken2_files)
      ch_versions    = ch_versions.mix(KRAKEN2.out.versions.first())
    }

    // TODO: Something clever with inital vs final qc

    INITIAL_QC(ch_initial_bam.map{ it -> tuple(it[0], it[1])})
    ch_versions    = ch_versions.mix(INITIAL_QC.out.versions.first())
    //ch_summary   = ch_summary.
    ch_for_multiqc = ch_for_multiqc.mix(INITIAL_QC.out.stats)

    FINAL_QC(ch_trim_bam.map{ it -> tuple(it[0], it[1])})
    ch_versions    = ch_versions.mix(FINAL_QC.out.versions.first())
    ch_for_multiqc = ch_for_multiqc.mix(FINAL_QC.out.flagstat)

    FINAL_QC.out.coverage
      .collectFile(name: "final_samtools_coverage_summary.tsv",
        keepHeader: true,
        storeDir: "${params.outdir}/samtools")
      .set { samtools_coverage_file }

    ch_summary = ch_summary.mix(FINAL_QC.out.stats)
    ch_summary = ch_summary.mix(FINAL_QC.out.flagstat)
    ch_summary = ch_summary.mix(samtools_coverage_file)
    ch_summary = ch_summary.mix(FINAL_QC.out.depth)

    if (params.aci) {
      ACI(ch_trim_bam.map{ it -> tuple(it[0], it[1])}.combine(ch_amplicon_bed))
      ch_versions    = ch_versions.mix(ACI.out.versions.first())
      ch_for_multiqc = ch_for_multiqc.mix(ACI.out.for_multiqc)
      
      ACI.out.cov
        .collectFile(name: "aci_coverage_summary.csv",
          keepHeader: true,
          storeDir: "${params.outdir}/aci")
        .set { aci_coverage_file }
      ch_summary = ch_summary.mix(aci_coverage_file)
    }

    if (params.bcftools_variants) {
      BCFTOOLS(ch_trim_bam.map{ it -> tuple(it[0], it[1])}.combine(ch_reference_genome))
      ch_versions = ch_versions.mix(BCFTOOLS.out.versions.first())
      ch_summary = ch_summary.mix(BCFTOOLS.out.bcftools_variants_file)

      if (params.igv_reports) {
        IGV_REPORTS(
          BCFTOOLS.out.vcf
          .join(ch_trim_bam, by: 0)
          .combine(ch_reference_genome)
        )
      ch_versions = ch_versions.mix(IGV_REPORTS.out.versions.first())
      }
    }

    if (params.ivar_variants) {

      IVAR_VARIANTS(ch_trim_bam.map{ it -> tuple(it[0], it[1])}.combine(ch_reference_genome).combine(ch_gff_file))
      ch_versions = ch_versions.mix(IVAR_VARIANTS.out.versions.first())
      ch_summary  = ch_summary.mix(IVAR_VARIANTS.out.variant_tsv)
    }

    if (params.samtools_ampliconstats) {
      AMPLICONSTATS(ch_trim_bam.map{ it -> tuple(it[0], it[1])}.combine(ch_primer_bed))
      ch_versions = ch_versions.mix(AMPLICONSTATS.out.versions.first())
      ch_summary  = ch_summary.mix(AMPLICONSTATS.out.ampliconstats.map{it -> it[1]})

      PLOT_AMPLICONSTATS(AMPLICONSTATS.out.ampliconstats)
      ch_versions = ch_versions.mix(PLOT_AMPLICONSTATS.out.versions.first())
    }



  emit:
    for_multiqc = ch_for_multiqc
    for_summary = ch_summary
    versions    = ch_versions
}
