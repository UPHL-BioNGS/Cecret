include { aci }                                                 from '../modules/aci'       addParams(params)
include { bcftools_variants }                                   from '../modules/bcftools'  addParams(params)
include { igv_reports }                                         from '../modules/igvreports' addParams(params)
include { ivar_variants }                                       from '../modules/ivar'      addParams(params)
include { fastqc }                                              from '../modules/fastqc'    addParams(params)
include { kraken2 }                                             from '../modules/kraken2'   addParams(params)
include { samtools_stats as samtools_intial_stats }             from '../modules/samtools'  addParams(params)
include { samtools_stats }                                      from '../modules/samtools'  addParams(params)
include { samtools_depth }                                      from '../modules/samtools'  addParams(params)
include { samtools_coverage }                                   from '../modules/samtools'  addParams(params)
include { samtools_ampliconstats; samtools_plot_ampliconstats } from '../modules/samtools'  addParams(params)
include { samtools_flagstat }                                   from '../modules/samtools'  addParams(params)

workflow qc {
  take:
    ch_raw_reads
    ch_clean_reads
    ch_kraken2_db
    ch_sam
    ch_trim_bam
    ch_reference_genome
    ch_gff_file
    ch_amplicon_bed
    ch_primer_bed

  main:
    ch_for_multiqc = Channel.empty()
    ch_for_summary = Channel.empty()
    ch_versions    = Channel.empty()
    
    fastqc(ch_raw_reads)
    ch_for_summary = ch_for_summary.mix(fastqc.out.fastq_name.collectFile(name: "fastq_names.csv", keepHeader: true ))
    ch_versions    = ch_versions.mix(fastqc.out.versions.first())

    if ( ch_kraken2_db ) {
      kraken2(ch_clean_reads.combine(ch_kraken2_db))
      ch_for_multiqc = ch_for_multiqc.mix(kraken2.out.kraken2_files)
      ch_for_summary = ch_for_summary.mix(kraken2.out.kraken2_files)
      ch_versions    = ch_versions.mix(kraken2.out.versions.first())
    }

    samtools_intial_stats( ch_sam)
    aci(                   ch_trim_bam.map{ it -> tuple(it[0], it[1])}.combine(ch_amplicon_bed))
    samtools_flagstat(     ch_trim_bam.map{ it -> tuple(it[0], it[1])})
    samtools_depth(        ch_trim_bam.map{ it -> tuple(it[0], it[1])})
    samtools_coverage(     ch_trim_bam.map{ it -> tuple(it[0], it[1])})
    samtools_stats(        ch_trim_bam.map{ it -> tuple(it[0], it[1])})
    bcftools_variants(     ch_trim_bam.map{ it -> tuple(it[0], it[1])}.combine(ch_reference_genome))
    ivar_variants(         ch_trim_bam.map{ it -> tuple(it[0], it[1])}.combine(ch_reference_genome).combine(ch_gff_file))
    samtools_ampliconstats(ch_trim_bam.map{ it -> tuple(it[0], it[1])}.combine(ch_primer_bed))
    samtools_plot_ampliconstats(samtools_ampliconstats.out.samtools_ampliconstats_files)

    ch_versions = ch_versions.mix(samtools_intial_stats.out.versions.first())
    ch_versions = ch_versions.mix(aci.out.versions.first())
    ch_versions = ch_versions.mix(samtools_flagstat.out.versions.first())
    ch_versions = ch_versions.mix(samtools_depth.out.versions.first())
    ch_versions = ch_versions.mix(samtools_coverage.out.versions.first())
    ch_versions = ch_versions.mix(samtools_stats.out.versions.first())
    ch_versions = ch_versions.mix(bcftools_variants.out.versions.first())
    ch_versions = ch_versions.mix(ivar_variants.out.versions.first())
    ch_versions = ch_versions.mix(samtools_ampliconstats.out.versions.first())
    ch_versions = ch_versions.mix(samtools_plot_ampliconstats.out.versions.first())

    bcftools_variants.out.vcf
      .join(ch_trim_bam, by: 0)
      .combine(ch_reference_genome)
      .set{ for_igv_reports }

    igv_reports(for_igv_reports)

    ch_versions = ch_versions.mix(igv_reports.out.versions.first())

    samtools_coverage.out.samtools_coverage
      .collectFile(name: "samtools_coverage_summary.tsv",
        keepHeader: true,
        storeDir: "${params.outdir}/samtools_coverage")
      .set { samtools_coverage_file }

    aci.out.cov
      .collectFile(name: "aci_coverage_summary.csv",
        keepHeader: true,
        storeDir: "${params.outdir}/aci")
      .set { aci_coverage_file }

    //# All of these tools are for QC, so each needs to make it to the summary file if actually useful
    ch_for_summary = ch_for_summary
      .mix(aci_coverage_file)
      .mix(ivar_variants.out.variant_tsv)
      .mix(bcftools_variants.out.bcftools_variants_file)
      .mix(samtools_stats.out.samtools_stats_files)
      .mix(samtools_depth.out.file)
      .mix(samtools_ampliconstats.out.samtools_ampliconstats_files.map{it -> it[1]})
      .mix(samtools_coverage_file)

  emit:
    for_multiqc = ch_for_multiqc.mix(fastqc.out.fastqc_files).mix(samtools_intial_stats.out.samtools_stats_files).mix(samtools_flagstat.out.samtools_flagstat_files).mix(aci.out.for_multiqc)
    for_summary = ch_for_summary
    versions    = ch_versions
}
