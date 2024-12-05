include { ARTIC                                 } from '../../modules/local/artic'
include { ARTIC_FILTER                          } from '../../modules/local/artic'
include { BBNORM                                } from '../../modules/local/bbnorm'
include { BWA                                   } from '../../modules/local/bwa'
include { FASTP                                 } from '../../modules/local/fastp'
include { IVAR_TRIM                             } from '../../modules/local/ivar'
include { IVAR_CONSENSUS as IVAR                } from '../../modules/local/ivar'
include { MINIMAP2                              } from '../../modules/local/minimap2'
include { SAMTOOLS_SORT as SORT                 } from '../../modules/local/samtools'
include { SAMTOOLS_AMPLICONCLIP as AMPLICONCLIP } from '../../modules/local/samtools'
include { SAMTOOLS_FILTER as FILTER             } from '../../modules/local/samtools'
include { SAMTOOLS_MARKDUP as MARKDUP           } from '../../modules/local/samtools'
include { SEQYCLEAN                             } from '../../modules/local/seqyclean'

workflow CONSENSUS {
  take:
    ch_reads
    ch_nanopore
    ch_reference
    ch_primer_bed

  main:
    ch_for_multiqc = Channel.empty()
    ch_for_version = Channel.empty()
    ch_versions    = Channel.empty()

  if ( params.bbnorm ){
    BBNORM(ch_reads)
    ch_norm_reads = BBNORM.out.fastq
    ch_versions   = ch_versions.mix(BBNORM.out.versions.first())
  } else {
    ch_norm_reads = ch_reads
  }

  if ( params.cleaner == 'seqyclean' ) {
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
    ch_for_version = ch_for_version.mix(SEQYCLEAN.out.cleaner_version)
    ch_for_multiqc = ch_for_multiqc.mix(SEQYCLEAN.out.seqyclean_files_collect_paired).mix(SEQYCLEAN.out.seqyclean_files_collect_single)
    ch_versions    = ch_versions.mix(SEQYCLEAN.out.versions.first())

  } else if ( params.cleaner == 'fastp' ) {
    FASTP(ch_norm_reads)

    ch_clean_reads = FASTP.out.clean_reads
    ch_for_version = ch_for_version.mix(FASTP.out.cleaner_version)
    ch_for_multiqc = ch_for_multiqc.mix(FASTP.out.fastp_files) 
    ch_versions    = ch_versions.mix(FASTP.out.versions.first())
  }

  if ( params.aligner == 'bwa' ) {
    BWA(ch_clean_reads.combine(ch_reference))

    ch_sam         = BWA.out.sam
    ch_for_version = ch_for_version.mix(BWA.out.aligner_version)
    ch_versions    = ch_versions.mix(BWA.out.versions.first())
  
  } else if ( params.aligner = 'minimap2' ) {
    MINIMAP2(ch_clean_reads.combine(ch_reference))
    
    ch_sam         = MINIMAP2.out.sam
    ch_for_version = ch_for_version.mix(MINIMAP2.out.aligner_version)
    ch_versions    = ch_versions.mix(MINIMAP2.out.versions.first())
  
  }

  if ( params.markdup ) {
    MARKDUP(ch_reads.join(ch_sam).map { it -> tuple(it[0], it[2], it[3])} )
    ch_bam = MARKDUP.out.bam_bai
    ch_versions = ch_versions.mix(MARKDUP.out.versions.first())
  
  } else {
    SORT(ch_sam)
    ch_bam = SORT.out.bam_bai
  
  }

  if ( params.trimmer == 'ivar' ) {
    IVAR_TRIM(ch_bam.filter{ it -> it[1].size() > 500 }.map{it -> tuple( it[0], it[1])}.combine(ch_primer_bed))

    ch_trim_bam    = IVAR_TRIM.out.bam_bai
    ch_for_version = ch_for_version.mix(IVAR_TRIM.out.trimmer_version)
    ch_for_multiqc = ch_for_multiqc.mix(IVAR_TRIM.out.ivar_trim_files)
    ch_versions    = ch_versions.mix(IVAR_TRIM.out.versions.first())

  } else if ( params.trimmer == 'samtools' || params.trimmer == 'ampliconclip' ) {
    AMPLICONCLIP(ch_bam.filter{ it -> it[1].size() > 500 }.map{it -> tuple( it[0], it[1])}.combine(ch_primer_bed))
    
    ch_trim_bam    = AMPLICONCLIP.out.bam_bai
    ch_for_version = ch_for_version.mix(AMPLICONCLIP.out.trimmer_version)
    ch_versions    = ch_versions.mix(AMPLICONCLIP.out.versions.first())
  
  } else if ( params.trimmer == 'none' ) {
    ch_trim_bam    = ch_bam
  
  }

  IVAR(ch_trim_bam.map{ it -> tuple(it[0], it[1])}.combine(ch_reference))
  ch_versions = ch_versions.mix(IVAR.out.versions.first())
  ch_for_version = ch_for_version.mix(IVAR.out.ivar_version)

  ARTIC_FILTER(ch_nanopore)
  ARTIC(ARTIC_FILTER.out.fastq.combine(ch_reference).combine(ch_primer_bed))
  ch_for_version = ch_for_version.mix(ARTIC.out.artic_version)
  ch_versions    = ch_versions.mix(ARTIC.out.versions.first()).mix(ARTIC_FILTER.out.versions.first())

  if ( params.filter ) { 
    FILTER(ch_sam.mix(ARTIC.out.bam))
    ch_filtered_reads = FILTER.out.filtered_reads
  } else {
    ch_filtered_reads = Channel.empty()
  }

  emit:
    consensus        = IVAR.out.consensus.mix(ARTIC.out.consensus)
    trim_bam         = ch_trim_bam.mix(ARTIC.out.bam)
    clean_reads      = ch_clean_reads
    sam              = ch_sam
    filtered_reads   = ch_filtered_reads

    for_multiqc      = ch_for_multiqc 
    for_version      = ch_for_version.unique().collect()

    versions         = ch_versions
}
