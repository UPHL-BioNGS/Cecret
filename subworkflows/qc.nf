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
      ch_for_summary = ch_for_summary.mix(kraken2.out.kraken2_files)
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

    ch_for_summary
      .mix(fastqc.out.fastqc_1_results)
      .mix(fastqc.out.fastqc_2_results)
      .mix(ivar_variants.out.ivar_variants_results)
      .mix(bcftools_variants.out.bcftools_variants_results)
      .mix(samtools_stats.out.samtools_stats_after_size_results)
      .mix(samtools_coverage.out.samtools_coverage_results)
      .mix(samtools_coverage.out.samtools_covdepth_results)
      .mix(samtools_depth.out.samtools_depth_results)
      .mix(samtools_ampliconstats.out.samtools_ampliconstats_results)
      .mix(bedtools_multicov.out.bedtools_results)
      .set (ch_all_summary)



  emit:
  for_multiqc = ch_for_multiqc.mix(fastqc.out.fastqc_files).mix(samtools_intial_stats.out.samtools_stats_files).mix(samtools_flagstat.out.samtools_flagstat_files)
  for_summary = ch_all_summary
}
