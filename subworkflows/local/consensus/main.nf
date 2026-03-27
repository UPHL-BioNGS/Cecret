include { ARTIC                                 } from '../../../modules/local/artic'
include { ARTIC_FILTER                          } from '../../../modules/local/artic_filter'
include { ARTIC_TOOLS                           } from '../../../modules/local/artic_tools'
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

  log.info """

Running consensus creation analysis. This workflow performs reference-based 
consensus generation and variant calling.

Relevant params and their values:
- 'params.primer_set' : ${params.primer_set}
    - Master variable with corresponding a corresponding primer and amplicon bed, 
      reference, and gff file.
- 'params.primer_bed' : ${params.primer_bed}
    - File designating the primers used to create the amplicons.
    - See https://github.com/UPHL-BioNGS/Cecret/wiki/Creating on creating a custom
      primer scheme.
- 'params.reference_genome' : ${params.reference_genome}
    - FASTA file used to align reads to.
    - Must match reference line in primer BED file.
    - See https://github.com/UPHL-BioNGS/Cecret/wiki/Reference for more information.
- 'params.cleaner' : ${params.cleaner}
    - Designates which tool will filter FASTQ files.
    - Options are 'seqyclean' and 'fastp'.
- 'params.aligner' : ${params.aligner}
    - Designates which tool will align clean reads to the reference.
    - OPtions are 'bwa' and 'minimap2'.
- 'params.trimmer' : ${params.trimmer}
    - Designates which tool will trim primers from aligned sequences.
    - Options are 'samtools', 'ivar', and 'none'.
    - Setting this as 'none' will skip primer trimming.
- 'params.filter' : ${params.filter}
    - Produces FASTQ files aligned to reference.
- 'params.markdup' : ${params.markdup}
    - Removes duplicates (may increase artifact rates in amplicon sequencing)

┏━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ process            ┃ description                                                       ┃
┣━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
┃ BBNORM             ┃ Normalizes large read datasets to a target depth.                 ┃
┃ SEQYCLEAN / FASTP  ┃ Filters out low-quality reads and removes adapters.               ┃
┃ BWA / MINIMAP2     ┃ Aligns the cleaned reads to the provided reference genome.        ┃
┃ SAMTOOLS           ┃ Sorts alignments, marks duplicates, and filters out non-targets.  ┃
┃ IVAR / AMPLICONCLIP┃ Trims amplicon primer sequences from the BAM alignments.          ┃
┃ IVAR CONSENSUS     ┃ Calls variants and generates the final consensus sequence.        ┃
┃ ARTIC              ┃ Filters Nanopore reads and generates Nanopore consensus fastas.   ┃
┃ ARTIC-TOOLS        ┃ Validates input primer bedfile for ARTIC.                         ┃
┗━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

"""


  ch_multiqc  = channel.empty()
  ch_versions = channel.empty()
  ch_bam      = channel.empty()
  ch_trim_bam = channel.empty()
  ch_nanopore_bam = channel.empty()
  ch_consensus = channel.empty()
  ch_clean_reads = channel.empty()

  // only show the following if there are illumina fastq files present
  if (params.sample_sheet || params.reads || ! params.sra_accessions.isEmpty()) {
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
    if ( params.cleaner == 'seqyclean') {
      // running seqyclean
      SEQYCLEAN(ch_norm_reads)

      SEQYCLEAN.out.seqyclean_files_collect_paired
        .collectFile(name: "Combined_SummaryStatistics.tsv",
          keepHeader: true,
          storeDir: "${params.outdir}/seqyclean")

      SEQYCLEAN.out.seqyclean_files_collect_single
        .collectFile(name: "Combined_seqyclean_SummaryStatistics.tsv",
          keepHeader: true,
          storeDir: "${params.outdir}/seqyclean")

      ch_clean_reads = ch_clean_reads.mix(SEQYCLEAN.out.clean_reads)
      ch_multiqc     = ch_multiqc.mix(SEQYCLEAN.out.seqyclean_files_collect_paired).mix(SEQYCLEAN.out.seqyclean_files_collect_single)
      ch_versions    = ch_versions.mix(SEQYCLEAN.out.versions.first())

    } else if ( params.cleaner == 'fastp' ) {
      // running fastp
      FASTP(ch_norm_reads)

      ch_clean_reads = ch_clean_reads.mix(FASTP.out.clean_reads)
      ch_multiqc     = ch_multiqc.mix(FASTP.out.fastp_files) 
      ch_versions    = ch_versions.mix(FASTP.out.versions.first())
    }

    // aligning reads to reference
    if ( params.aligner == 'bwa' ) {
      // running bwa
      BWA(ch_clean_reads.combine(ch_reference))

      ch_sam      = BWA.out.sam
      ch_versions = ch_versions.mix(BWA.out.versions.first())
    
    } else if ( params.aligner == 'minimap2') {
      // running minimap2
      MINIMAP2(ch_clean_reads.combine(ch_reference))
      
      ch_sam      = MINIMAP2.out.sam
      ch_versions = ch_versions.mix(MINIMAP2.out.versions.first())
    
    } else {
      ch_sam = channel.empty()
    }

    // removing duplicates
    if ( params.markdup ) {
      MARKDUP(ch_reads.join(ch_sam).map { it -> tuple(it[0], it[2], it[3])} )
      ch_bam      = ch_bam.mix(MARKDUP.out.bam_bai)
      ch_versions = ch_versions.mix(MARKDUP.out.versions.first())
    
    } else {
      SORT(ch_sam)
      ch_bam = ch_bam.mix(SORT.out.bam_bai)
      ch_versions = ch_versions.mix(SORT.out.versions.first())
    }

    // removing primers
    if ( params.trimmer == 'ivar' ) {
      // running ivar trim
      IVAR_TRIM(ch_bam.filter{ it -> it[1].size() > 500 }.map{it -> tuple( it[0], it[1])}.combine(ch_primer_bed))

      ch_trim_bam = ch_trim_bam.mix(IVAR_TRIM.out.bam_bai)
      ch_multiqc  = ch_multiqc.mix(IVAR_TRIM.out.ivar_trim_files)
      ch_versions = ch_versions.mix(IVAR_TRIM.out.versions.first())

    } else if ( params.trimmer == 'samtools' || params.trimmer == 'ampliconclip' ) {
      // running samtools ampliconclip
      AMPLICONCLIP(ch_bam.filter{ it -> it[1].size() > 500 }.map{it -> tuple( it[0], it[1])}.combine(ch_primer_bed))
      
      ch_trim_bam = ch_trim_bam.mix(AMPLICONCLIP.out.bam_bai)
      ch_versions = ch_versions.mix(AMPLICONCLIP.out.versions.first())
    
    } else if ( params.trimmer == 'none' ) {
      // mostly for bait-derived libraries
      // skipping trimming (not recommended for amplicon-created libraries)
      ch_trim_bam = ch_trim_bam.mix(ch_bam)
    }

    // getting a consensus with ivar
    IVAR(ch_trim_bam.map{ it -> tuple(it[0], it[1])}.combine(ch_reference))
    ch_consensus = ch_consensus.mix(IVAR.out.consensus)
    ch_versions = ch_versions.mix(IVAR.out.versions.first())
  }
  // running artic on nanopore reads
  // only show if there is a sample sheet as input (which could have nanopor reads),
  // or nanopore reads are read in through a channel
  if (params.sample_sheet || params.nanopore) {
    if (params.artic && params.artic_filter) {
      ARTIC_FILTER(ch_nanopore)
      ch_versions = ch_versions.mix(ARTIC_FILTER.out.versions.first())

      if (params.primer_bed) {
        ARTIC_TOOLS(ch_primer_bed)
        ch_versions = ch_versions.mix(ARTIC_TOOLS.out.versions.first())
        ch_primer = ARTIC_TOOLS.out.bed
      } else {
        ch_primer = ch_primer_bed
      }

      if (params.artic) {
        ARTIC(ARTIC_FILTER.out.fastq.combine(ch_reference).combine(ch_primer))
        ch_versions = ch_versions.mix(ARTIC.out.versions.first())
        ch_nanopore_bam = ch_nanopore_bam.mix(ARTIC.out.bam)
        ch_consensus = ch_consensus.mix(ARTIC.out.consensus)
      }
    }
  }

  // remove all non-target-organism reads
  // simpler for non-human hosts or when more than one organism needs to be removed
  if ( params.filter ) { 
    FILTER(ch_sam.mix(ch_nanopore_bam))
    ch_filtered_reads = FILTER.out.filtered_reads
    ch_versions       = ch_versions.mix(FILTER.out.versions.first())
  } else {
    ch_filtered_reads = channel.empty()
  }

  emit:
    consensus        = ch_consensus
    just_bam         = ch_bam
    trim_bam         = ch_trim_bam.mix(ch_nanopore_bam)
    clean_reads      = ch_clean_reads
    filtered_reads   = ch_filtered_reads
    for_multiqc      = ch_multiqc 
    versions         = ch_versions
}
