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
}
