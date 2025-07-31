include { ARTIC                                 } from '../../../modules/local/artic'
include { ARTIC_FILTER                          } from '../../../modules/local/artic'
include { BBNORM                                } from '../../../modules/local/bbnorm'
include { BWA                                   } from '../../../modules/local/bwa/'
include { FASTP                                 } from '../../../modules/local/fastp'
include { IVAR_TRIM                             } from '../../../modules/local/ivar'
include { IVAR_CONSENSUS as IVAR                } from '../../../modules/local/ivar'
include { MINIMAP2                              } from '../../../modules/local/minimap2'
include { SAMTOOLS_SORT as SORT                 } from '../../../modules/local/samtools'
include { SAMTOOLS_AMPLICONCLIP as AMPLICONCLIP } from '../../../modules/local/samtools'
include { SAMTOOLS_FILTER as FILTER             } from '../../../modules/local/samtools'
include { SAMTOOLS_MARKDUP as MARKDUP           } from '../../../modules/local/samtools'
include { SEQYCLEAN                             } from '../../../modules/local/seqyclean'

workflow CONSENSUS {
  take:
  ch_reads // channel: [meta, reads]
  ch_nanopore // channel: [meta, reads]
  ch_reference // channel: fasta
  ch_primer_bed // channel: bedfile

  main:
  ch_multiqc  = Channel.empty()
  ch_versions = Channel.empty()

  // running bbnorm to normalize large datasets
  // this is not recommended for wastewater
  if ( params.bbnorm ){
    BBNORM(ch_reads)
    ch_norm_reads = BBNORM.out.fastq
    ch_versions   = ch_versions.mix(BBNORM.out.versions.first())
  } else {
    ch_norm_reads = ch_reads
  }

  // filtering out low quality reads
  if ( params.cleaner == 'seqyclean' ) {
    // running seqyclean
    SEQYCLEAN(ch_norm_reads)

    SEQYCLEAN.out.seqyclean_files_collect_paired
      .collectFile(name: "Combined_SummaryStatistics.tsv",
        keepHeader: true,
        storeDir: "${params.outdir}/seqyclean")
      .set { seqyclean_file1 }

    SEQYCLEAN.out.seqyclean_files_collect_single
      .collectFile(name: "Combined_seqyclean_SummaryStatistics.tsv",
        keepHeader: true,
        storeDir: "${params.outdir}/seqyclean")
      .set { seqyclean_file2 }

    ch_clean_reads = SEQYCLEAN.out.clean_reads
    ch_multiqc     = ch_multiqc.mix(SEQYCLEAN.out.seqyclean_files_collect_paired).mix(SEQYCLEAN.out.seqyclean_files_collect_single)
    ch_versions    = ch_versions.mix(SEQYCLEAN.out.versions.first())

  } else if ( params.cleaner == 'fastp' ) {
    // running fastp
    FASTP(ch_norm_reads)

    ch_clean_reads = FASTP.out.clean_reads
    ch_multiqc     = ch_multiqc.mix(FASTP.out.fastp_files) 
    ch_versions    = ch_versions.mix(FASTP.out.versions.first())
  }

  // aligning reads to reference
  if ( params.aligner == 'bwa' ) {
    // running bwa
    BWA(ch_clean_reads.combine(ch_reference))

    ch_sam      = BWA.out.sam
    ch_versions = ch_versions.mix(BWA.out.versions.first())
  
  } else if ( params.aligner == 'minimap2' ) {
    // running minimap2
    MINIMAP2(ch_clean_reads.combine(ch_reference))
    
    ch_sam      = MINIMAP2.out.sam
    ch_versions = ch_versions.mix(MINIMAP2.out.versions.first())
  
  }

  // removing duplicates
  if ( params.markdup ) {
    MARKDUP(ch_reads.join(ch_sam).map { it -> tuple(it[0], it[2], it[3])} )
    ch_bam      = MARKDUP.out.bam_bai
    ch_versions = ch_versions.mix(MARKDUP.out.versions.first())
  
  } else {
    SORT(ch_sam)
    ch_bam = SORT.out.bam_bai
    ch_versions = ch_versions.mix(SORT.out.versions.first())
  }

  // removing primers
  if ( params.trimmer == 'ivar' ) {
    // running ivar trim
    IVAR_TRIM(ch_bam.filter{ it -> it[1].size() > 500 }.map{it -> tuple( it[0], it[1])}.combine(ch_primer_bed))

    ch_trim_bam = IVAR_TRIM.out.bam_bai
    ch_multiqc  = ch_multiqc.mix(IVAR_TRIM.out.ivar_trim_files)
    ch_versions = ch_versions.mix(IVAR_TRIM.out.versions.first())

  } else if ( params.trimmer == 'samtools' || params.trimmer == 'ampliconclip' ) {
    // running samtools ampliconclip
    AMPLICONCLIP(ch_bam.filter{ it -> it[1].size() > 500 }.map{it -> tuple( it[0], it[1])}.combine(ch_primer_bed))
    
    ch_trim_bam = AMPLICONCLIP.out.bam_bai
    ch_versions = ch_versions.mix(AMPLICONCLIP.out.versions.first())
  
  } else if ( params.trimmer == 'none' ) {
    // mostly for bait-derived libraries
    // skipping trimming (not recommended for amplicon-created libraries)
    ch_trim_bam = ch_bam
  }

  // getting a consensus with ivar
  IVAR(ch_trim_bam.map{ it -> tuple(it[0], it[1])}.combine(ch_reference))
  ch_versions = ch_versions.mix(IVAR.out.versions.first())

  // running artic on nanopore reads
  // TODO : Hide this if there are no nanopore reads
  if (params.artic && params.artic_filter) {
    ARTIC_FILTER(ch_nanopore)
    ch_versions = ch_versions.mix(ARTIC_FILTER.out.versions.first())

    if (params.artic) {
      ARTIC(ARTIC_FILTER.out.fastq.combine(ch_reference).combine(ch_primer_bed))
      ch_versions = ch_versions.mix(ARTIC.out.versions.first())
    }
  }

  // remove all non-target-organism reads
  // simpler for non-human hosts or when more than one organism needs to be removed
  if ( params.filter ) { 
    FILTER(ch_sam.mix(ARTIC.out.bam))
    ch_filtered_reads = FILTER.out.filtered_reads
    ch_versions       = ch_versions.mix(FILTER.out.versions.first())
  } else {
    ch_filtered_reads = Channel.empty()
  }

  emit:
    consensus        = IVAR.out.consensus.mix(ARTIC.out.consensus)
    just_bam         = ch_bam
    trim_bam         = ch_trim_bam.mix(ARTIC.out.bam)
    clean_reads      = ch_clean_reads
    filtered_reads   = ch_filtered_reads

    for_multiqc      = ch_multiqc 

    versions         = ch_versions
}
