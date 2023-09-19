include { artic ; artic_read_filtering }          from '../modules/artic'     addParams(params)
include { bwa }                                   from '../modules/bwa'       addParams(params)
include { fastp }                                 from '../modules/fastp'     addParams(params)
include { ivar_trim ; ivar_consensus as ivar }    from '../modules/ivar'      addParams(params)
include { minimap2 }                              from '../modules/minimap2'  addParams(params)
include { samtools_sort as sort }                 from '../modules/samtools'  addParams(params)
include { samtools_ampliconclip as ampliconclip } from '../modules/samtools'  addParams(params)
include { samtools_filter as filter }             from '../modules/samtools'  addParams(params)
include { samtools_markdup as markdup }           from '../modules/samtools'  addParams(params)
include { seqyclean }                             from '../modules/seqyclean' addParams(params)

workflow cecret {
  take:
    ch_reads
    ch_nanopore
    ch_reference
    ch_primer_bed

  main:
    ch_for_multiqc = Channel.empty()
    ch_for_version = Channel.empty()

  if ( params.cleaner == 'seqyclean' ) {
    seqyclean(ch_reads)

    seqyclean.out.seqyclean_files_collect_paired
      .collectFile(name: "Combined_SummaryStatistics.tsv",
        keepHeader: true,
        storeDir: "${params.outdir}/seqyclean")
      .set { seqyclean_file1 }

    seqyclean.out.seqyclean_files_collect_single
      .collectFile(name: "Combined_seqyclean_SummaryStatistics.tsv",
        keepHeader: true,
        storeDir: "${params.outdir}/seqyclean")
      .set { seqyclean_file2 }

    ch_clean_reads = seqyclean.out.clean_reads
    ch_for_version = ch_for_version.mix(seqyclean.out.cleaner_version)
    ch_for_multiqc = ch_for_multiqc.mix(seqyclean.out.seqyclean_files_collect_paired).mix(seqyclean.out.seqyclean_files_collect_single)

  } else if ( params.cleaner == 'fastp' ) {
    fastp(ch_reads)

    ch_clean_reads = fastp.out.clean_reads
    ch_for_version = ch_for_version.mix(fastp.out.cleaner_version)
    ch_for_multiqc = ch_for_multiqc.mix(fastp.out.fastp_files) 
  }

  if ( params.aligner == 'bwa' ) {
    bwa(ch_clean_reads.combine(ch_reference))

    ch_sam         = bwa.out.sam
    ch_for_version = ch_for_version.mix(bwa.out.aligner_version)
  
  } else if ( params.aligner = 'minimap2' ) {
    minimap2(ch_clean_reads.combine(ch_reference))
    
    ch_sam         = minimap2.out.sam
    ch_for_version = ch_for_version.mix(minimap2.out.aligner_version)
  
  }

  if ( params.markdup ) {
    markdup(ch_reads.join(ch_sam).map { it -> tuple(it[0], it[2], it[3])} )
    ch_bam = markdup.out.bam_bai
  
  } else {
    sort(ch_sam)
    ch_bam = sort.out.bam_bai
  
  }

  if ( params.trimmer == 'ivar' ) {
    ivar_trim(ch_bam.filter{ it -> it[1].size() > 500 }.map{it -> tuple( it[0], it[1])}.combine(ch_primer_bed))

    ch_trim_bam    = ivar_trim.out.bam_bai
    ch_for_version = ch_for_version.mix(ivar_trim.out.trimmer_version)
    ch_for_multiqc = ch_for_multiqc.mix(ivar_trim.out.ivar_trim_files)

  } else if ( params.trimmer == 'samtools' || params.trimmer == 'ampliconclip' ) {
    ampliconclip(ch_bam.filter{ it -> it[1].size() > 500 }.map{it -> tuple( it[0], it[1])}.combine(ch_primer_bed))
    
    ch_trim_bam    = ampliconclip.out.bam_bai
    ch_for_version = ch_for_version.mix(ampliconclip.out.trimmer_version)
  
  } else if ( params.trimmer == 'none' ) {
    ch_trim_bam    = ch_bam
  
  }

  ivar(ch_trim_bam.map{ it -> tuple(it[0], it[1])}.combine(ch_reference))

  artic_read_filtering(ch_nanopore)
  artic(artic_read_filtering.out.fastq.combine(ch_reference).combine(ch_primer_bed))
  ch_for_version = ch_for_version.mix(artic.out.artic_version).mix(ivar.out.ivar_version)

  if ( params.filter ) { filter(ch_sam.mix(artic.out.bam)) }

  emit:
    consensus        = ivar.out.consensus.mix(artic.out.consensus)
    trim_bam         = ch_trim_bam.mix(artic.out.bam)
    clean_reads      = ch_clean_reads
    sam              = ch_sam

    for_multiqc      = ch_for_multiqc 
    for_version      = ch_for_version.unique().collect()
}
