include { bcftools_variants }                                   from '../modules/bcftools'  addParams(bcftools_variants_options: params.bcftools_variants_options )
include { ivar_variants }                                       from '../modules/ivar'      addParams(ivar_variants_options: params.ivar_variants_options)
include { bedtools_multicov }                                   from '../modules/bedtools'  addParams(bedtools_multicov_options: params.bedtools_multicov_options)
include { fastqc }                                              from '../modules/fastqc'    addParams(fastqc_options: params.fastqc_options )
include { kraken2 }                                             from '../modules/kraken2'   addParams(kraken2_options: params.kraken2_options, kraken2_organism: params.kraken2_organism)
include { samtools_stats as samtools_intial_stats }             from '../modules/samtools'  addParams(samtools_stats: params.samtools_stats, samtools_stats_options: params.samtools_stats_options)
include { samtools_stats }                                      from '../modules/samtools'  addParams(samtools_stats: params.samtools_stats, samtools_stats_options: params.samtools_stats_options)
include { samtools_depth }                                      from '../modules/samtools'  addParams(samtools_depth: params.samtools_depth, samtools_depth_options: params.samtools_depth_options)
include { samtools_coverage }                                   from '../modules/samtools'  addParams(samtools_coverage: params.samtools_coverage, samtools_coverage_options: params.samtools_coverage_options)
include { samtools_ampliconstats; samtools_plot_ampliconstats } from '../modules/samtools'  addParams(samtools_ampliconstats: params.samtools_ampliconstats, samtools_ampliconstats_options: params.samtools_ampliconstats_options, samtools_plot_ampliconstats: params.samtools_plot_ampliconstats, samtools_plot_ampliconstats_options: params.samtools_plot_ampliconstats_options)
include { samtools_flagstat }                                   from '../modules/samtools'  addParams(samtools_flagstat: params.samtools_flagstat, samtools_flagstat_options: params.samtools_flagstat_options)
include { freyja; freyja_aggregate }                            from '../modules/freyja'    addParams(freyja: params.freyja, freyja_variants_options: params.freyja_variants_options, freyja_demix_options: params.freyja_demix_options, freyja_aggregate: params.freyja_aggregate, freyja_aggregate_options: params.freyja_aggregate_options, freyja_plot_options: params.freyja_plot_options)

workflow qc {
  take:
  reads
  clean_reads
  kraken2_db
  sam
  bam
  bam_bai
  reference_genome
  gff_file
  amplicon_bed
  primer_bed

  main:
  fastqc(reads)
  if ( kraken2_db ) {
    kraken2(clean_reads.combine(kraken2_db))
    kraken2_files = kraken2.out.kraken2_files
  } else {
    kraken2_files = Channel.empty()
  }
  samtools_flagstat(bam)
  samtools_depth(bam)
  samtools_coverage(bam)
  samtools_stats(bam)
  samtools_intial_stats(sam)
  samtools_ampliconstats(bam.combine(primer_bed))

  samtools_plot_ampliconstats(samtools_ampliconstats.out.samtools_ampliconstats_files)

  bcftools_variants(bam.combine(reference_genome))
  ivar_variants(bam.combine(reference_genome).combine(gff_file))
  bedtools_multicov(bam_bai.combine(amplicon_bed))

  freyja(bam.combine(reference_genome))
  freyja_aggregate(freyja.out.freyja_demix.collect())

  emit:
  // for multiqc
  fastqc_files                    = fastqc.out.fastqc_files
  samtools_stats_files            = samtools_intial_stats.out.samtools_stats_files
  kraken2_files                   = kraken2_files
  samtools_flagstat_files         = samtools_flagstat.out.samtools_flagstat_files

  // for summary file
  fastqc_1_results                = fastqc.out.fastqc_1_results
  fastqc_2_results                = fastqc.out.fastqc_2_results
  kraken2_target_results          = kraken2.out.kraken2_target_results
  kraken2_human_results           = kraken2.out.kraken2_human_results
  ivar_variants_results           = ivar_variants.out.ivar_variants_results
  bcftools_variants_results       = bcftools_variants.out.bcftools_variants_results
  insert_size_after_trimming      = samtools_stats.out.samtools_stats_after_size_results
  samtools_coverage_results       = samtools_coverage.out.samtools_coverage_results
  samtools_covdepth_results       = samtools_coverage.out.samtools_covdepth_results
  samtools_depth_results          = samtools_depth.out.samtools_depth_results
  samtools_ampliconstats_results  = samtools_ampliconstats.out.samtools_ampliconstats_results
  bedtools_results                = bedtools_multicov.out.bedtools_results
  freyja_file                     = freyja_aggregate.out.aggregated_freyja_file
}
