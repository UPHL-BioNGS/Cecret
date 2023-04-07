include { bcftools_variants }                                   from '../modules/bcftools'  addParams(params)
include { ivar_variants }                                       from '../modules/ivar'      addParams(params)
include { bedtools_multicov }                                   from '../modules/bedtools'  addParams(params)
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
    
    fastqc(ch_raw_reads)
    if ( ch_kraken2_db ) {
      kraken2(ch_clean_reads.combine(ch_kraken2_db))
      ch_for_multiqc = ch_for_multiqc.mix(kraken2.out.kraken2_files)
    } 

    samtools_intial_stats( ch_sam)
    samtools_flagstat(     ch_trim_bam.map{ it -> tuple(it[0], it[1])})
    samtools_depth(        ch_trim_bam.map{ it -> tuple(it[0], it[1])})
    samtools_coverage(     ch_trim_bam.map{ it -> tuple(it[0], it[1])})
    samtools_stats(        ch_trim_bam.map{ it -> tuple(it[0], it[1])})
    bcftools_variants(     ch_trim_bam.map{ it -> tuple(it[0], it[1])}.combine(ch_reference_genome))
    ivar_variants(         ch_trim_bam.map{ it -> tuple(it[0], it[1])}.combine(ch_reference_genome).combine(ch_gff_file))
    samtools_ampliconstats(ch_trim_bam.map{ it -> tuple(it[0], it[1])}.combine(ch_primer_bed))
    samtools_plot_ampliconstats(samtools_ampliconstats.out.samtools_ampliconstats_files)
    bedtools_multicov(ch_trim_bam.combine(ch_amplicon_bed))

    samtools_coverage.out.samtools_coverage
      .collectFile(name: "samtools_coverage_summary.tsv",
        keepHeader: true,
        storeDir: "${params.outdir}/samtools_coverage")
      .set { samtools_coverage_file }


    //# All of these tools are for QC, so each needs to make it to the summary file if actually useful
    ch_for_summary = ch_for_summary
      .mix(ivar_variants.out.variant_tsv)
      .mix(bcftools_variants.out.bcftools_variants_file)
      .mix(samtools_stats.out.samtools_stats_files)
      .mix(samtools_depth.out.file)
      .mix(samtools_ampliconstats.out.samtools_ampliconstats_files)
      .mix(bedtools_multicov.out.multicov)

  emit:
    for_multiqc = ch_for_multiqc.mix(fastqc.out.fastqc_files).mix(samtools_intial_stats.out.samtools_stats_files).mix(samtools_flagstat.out.samtools_flagstat_files).mix(samtools_coverage_file)
    for_summary = ch_for_summary
}
