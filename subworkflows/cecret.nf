include { seqyclean }                                                                                                             from '../modules/seqyclean' addParams(params)
include { fastp }                                                                                                                 from '../modules/fastp'     addParams(params)
include { bwa }                                                                                                                   from '../modules/bwa'       addParams(params)
include { minimap2 }                                                                                                              from '../modules/minimap2'  addParams(params)
include { ivar_trim ; ivar_consensus as ivar }                                                                                    from '../modules/ivar'      addParams(params)
include { samtools_sort as sort; samtools_ampliconclip as ampliconclip; samtools_filter as filter; samtools_markdup as markdup }  from '../modules/samtools'  addParams(params)

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
    seqyclean_files1  = seqyclean.out.seqyclean_files_collect_paired
    seqyclean_files2  = seqyclean.out.seqyclean_files_collect_single
    clean_type        = seqyclean.out.clean_reads
    fastp_results     = Channel.empty()
    fastp_files       = Channel.empty()
  } else if ( params.cleaner == 'fastp' ) {
    fastp(reads)
    clean_reads       = fastp.out.paired_files.mix(fastp.out.single_files)
    clean_type        = fastp.out.clean_reads
    cleaner_version   = fastp.out.cleaner_version
    fastp_results     = fastp.out.fastp_results
    fastp_files       = fastp.out.fastp_files
    seqyclean_files1  = Channel.empty()
    seqyclean_files2  = Channel.empty()
  }

  if ( params.aligner == 'bwa' ) {
    bwa(clean_reads.combine(reference))
    sam               = bwa.out.sam
    aligner_version   = bwa.out.aligner_version
  } else if ( params.aligner = 'minimap2' ) {
    minimap2(clean_reads.combine(reference))
    sam               = minimap2.out.sam
    aligner_version   = minimap2.out.aligner_version
  }

  if ( params.markdup ) {
    markdup(reads.join(sam).map { reads -> tuple(reads[0], reads[2], reads[3])} )
    bam = markdup.out.bam
    bam_bai = markdup.out.bam_bai
  } else {
    sort(sam)
    bam = sort.out.bam
    bam_bai = sort.out.bam_bai
  }

  if ( params.trimmer == 'ivar' ) {
    ivar_trim(bam.combine(primer_bed))
    trimmed_bam       = ivar_trim.out.trimmed_bam
    trimmer_version   = ivar_trim.out.trimmer_version
    bam_bai           = ivar_trim.out.bam_bai
    ivar_files        = ivar_trim.out.ivar_trim_files
  } else if ( params.trimmer == 'samtools' ) {
    ampliconclip(bam.combine(primer_bed))
    trimmed_bam       = ampliconclip.out.trimmed_bam
    bam_bai           = ampliconclip.out.bam_bai
    trimmer_version   = ampliconclip.out.trimmer_version
    ivar_files        = Channel.empty()
  } else if ( params.trimmer == 'none' ) {
    trimmed_bam       = bam
    trimmer_version   = reads.map { reads -> tuple(reads[0],'none') }
    bam_bai           = bam_bai
    ivar_files        = Channel.empty()
  }

  ivar(trimmed_bam.combine(reference))

  if ( params.filter ) { filter(sam) }

  emit:
  consensus         = ivar.out.consensus
  bam               = trimmed_bam
  bam_bai           = bam_bai
  clean_type        = clean_type
  sam               = sam

  // for multiqc
  seqyclean_files1  = seqyclean_files1
  seqyclean_files2  = seqyclean_files2
  fastp_files       = fastp_files
  ivar_files        = ivar_files

  // for summary file
  consensus_results = ivar.out.consensus_results
  cleaner_version   = cleaner_version
  aligner_version   = aligner_version
  trimmer_version   = trimmer_version
  ivar_version      = ivar.out.ivar_version
  fastp_results     = fastp_results
}
