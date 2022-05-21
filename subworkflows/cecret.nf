include { seqyclean }                                                                                 from '../modules/seqyclean' addParams(seqyclean_options: params.seqyclean_options )
include { fastp }                                                                                     from '../modules/fastp'     addParams(fastp_options: params.fastp_options)
include { bwa }                                                                                       from '../modules/bwa'       addParams()
include { minimap2 }                                                                                  from '../modules/minimap2'  addParams(minimap2_options: params.minimap2_options )
include { ivar_trim ; ivar_consensus as ivar }                                                        from '../modules/ivar'      addParams(ivar_trim_options: params.ivar_trim_options)
include { samtools_sort as sort; samtools_ampliconclip as ampliconclip; samtools_filter as filter }   from '../modules/samtools'  addParams(samtools_ampliconclip_options: params.samtools_ampliconclip_options )

workflow cecret {
  take:
  reads
  reference
  primer_bed

  main:
  if ( params.cleaner == 'seqyclean' ) {
    seqyclean(reads)
    clean_reads       = seqyclean.out.single_reads.mix(seqyclean.out.paired_reads)
    cleaner_version   = seqyclean.out.cleaner_version
    seqyclean_files   = seqyclean.out.seqyclean_files_collect
    clean_type        = seqyclean.out.clean_reads
    fastp_results     = Channel.empty()
    fastp_files       = Channel.empty()
  } else if ( params.cleaner == 'fastp' ) {
    fastp(reads)
    clean_reads       = fastp.out.reads
    clean_type        = fastp.out.clean_type
    cleaner_version   = fastp.out.cleaner_version
    fastp_results     = fastp.out.fastp_results
    fastp_files       = fastp.out.fastp_files
    seqyclean_files   = Channel.empty()
  }

  if ( params.aligner == 'bwa' ) {
    bwa(clean_reads.combine(reference))
    sam               = bwa.out.sam
    aligner_version   = bwa.out.aligner_version
  } else if ( params.aligner = 'minimap2' ) {
    minimap2(clean_reads, reference)
    sam               = minimap2.out.sam
    aligner_version   = minimap2.out.aligner_version
  }
  sort(sam)

  if ( params.trimmer == 'ivar' ) {
    ivar_trim(sort.out.bam.combine(primer_bed))
    trimmed_bam       = ivar_trim.out.trimmed_bam
    trimmer_version   = ivar_trim.out.trimmer_version
    bam_bai           = ivar_trim.out.bam_bai
  } else if ( params.trimmer == 'samtools' ) {
    ampliconclip(sort.out.bam.combine(primer_bed))
    trimmed_bam       = ampliconclip.out.trimmed_bam
    bam_bai           = ampliconclip.out.bam_bai
    trimmer_version   = ampliconclip.out.trimmer_version
  } else if ( params.trimmer == 'none' ) {
    trimmed_bam       = sort.out.bam
    trimmer_version   = 'none'
    bam_bai           = sort.out.bam_bai
  }

  ivar(trimmed_bam.combine(reference))

  if ( params.filter ) { filter(sam) }

  emit:
  consensus         = ivar.out.consensus
  bam               = trimmed_bam
  bam_bai           = bam_bai
  clean_type        = clean_type

  // for multiqc
  seqyclean_files   = seqyclean_files
  fastp_files       = fastp_files
  ivar_files        = ivar_trim.out.ivar_trim_files

  // for summary file
  fastp_results     = fastp_results
  cleaner_version   = cleaner_version
  aligner_version   = aligner_version
  trimmer_version   = trimmer_version
  ivar_version      = ivar.out.ivar_version
  consensus_results = ivar.out.consensus_results
}
