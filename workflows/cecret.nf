/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SUMMARY       } from '../modules/local/local'
include { CONSENSUS     } from '../subworkflows/local/consensus'
include { QC            } from '../subworkflows/local/qc'
include { MSA           } from '../subworkflows/local/msa'
include { MULTIQC       } from '../modules/local/multiqc'
include { MPOX          } from '../subworkflows/local/mpox'                             
include { MPOX as OTHER } from '../subworkflows/local/mpox'
include { SARSCOV2      } from '../subworkflows/local/sarscov2'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CECRET {

    take:
    ch_reads // channel: [id, reads]
    ch_fastas // channel: [id, fasta]
    ch_multifasta // channel: fasta
    ch_nanopore // channel: [id, reads]
    ch_reference // channel: fasta
    ch_gff // channel: gff file
    ch_primer // channel: bedfile
    ch_amplicon // channel: bedfile
    ch_versions // channel: value
    ch_kraken2_db // channel: path
    ch_scripts // channel: [scripts]
    ch_nextclade_dataset // channel: file

    main:
    ch_for_version = Channel.empty()

    CONSENSUS(
        ch_reads,
        ch_nanopore,
        ch_reference,
        ch_primer
    )
    ch_versions = ch_versions.mix(CONSENSUS.out.versions)
    ch_consensus = CONSENSUS.out.consensus.mix(ch_fastas).mix(ch_multifasta)
    ch_for_version = CONSENSUS.out.for_version

    QC(ch_reads,
      CONSENSUS.out.clean_reads,
      ch_kraken2_db,
      CONSENSUS.out.sam,
      CONSENSUS.out.trim_bam,
      ch_reference,
      ch_gff,
      ch_amplicon,
      ch_primer)

    ch_for_multiqc = CONSENSUS.out.for_multiqc.mix(QC.out.for_multiqc)
    ch_for_summary = QC.out.for_summary
    ch_versions    = ch_versions.mix(QC.out.versions)

    if ( params.species == 'sarscov2' ) {
      SARSCOV2(
        ch_consensus,
        cecret.out.trim_bam,
        ch_reference,
        ch_nextclade_dataset,
        ch_scripts
        )

      ch_for_multiqc = ch_for_multiqc.mix(SARSCOV2.out.for_multiqc)
      ch_for_summary = ch_for_summary.mix(SARSCOV2.out.for_summary)
      ch_versions    = ch_versions.mix(SARSCOV2.out.versions)
    
    } else if ( params.species == 'mpx') {
      MPOX(ch_consensus)
      
      ch_for_multiqc = ch_for_multiqc.mix(MPOX.out.for_multiqc)
      ch_for_summary = ch_for_summary.mix(MPOX.out.for_summary)
      ch_versions    = ch_versions.mix(MPOX.out.versions)

    } else if ( params.species == 'other') {
      OTHER(ch_consensus)
      
      ch_for_multiqc = ch_for_multiqc.mix(OTHER.out.for_multiqc)
      ch_for_summary = ch_for_summary.mix(OTHER.out.for_summary)
      ch_versions    = ch_versions.mix(OTHER.out.versions)

    } 

    if ( params.relatedness ) { 
      
      MSA(
        ch_consensus,
        ch_reference
        ) 

      tree      = MSA.out.tree
      alignment = MSA.out.msa
      matrix    = MSA.out.matrix

      ch_for_multiqc = ch_for_multiqc.mix(MSA.out.for_multiqc)
      ch_versions    = ch_versions.mix(MSA.out.versions)

    } else {
      tree      = Channel.empty()
      alignment = Channel.empty()
      matrix    = Channel.empty()
    }

    ch_versions
      .collectFile(
        keepHeader: false,
        name: "collated_versions.yml")
      .set { ch_collated_versions }

    MULTIQC(
        ch_for_multiqc.mix(ch_collated_versions).collect(), 
        ch_scripts
    )
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)

    SUMMARY(
      ch_for_summary.mix(ch_fastas).mix(CONSENSUS.out.consensus).collect().map{it -> tuple([it])}
        .combine(ch_scripts)
        .combine(ch_for_version.mix(CONSENSUS.out.for_version).collect().map{it -> tuple([it])})
        .combine(MULTIQC.out.files.ifEmpty([]).map{it -> tuple([it])}))

  emit:
    bam       = CONSENSUS.out.trim_bam
    consensus = ch_consensus
    tree      = tree
    alignment = alignment
    matrix    = matrix
}