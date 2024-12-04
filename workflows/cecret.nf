/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SUMMARY  } from './modules/cecret'
include { LOCAL                } from './subworkflows/cecret'
include { QC }                   from './subworkflows/qc'
include { MSA }                  from './subworkflows/msa'
include { MULTIQC } from './modules/multiqc'
include { MPX }                  from './subworkflows/mpx'                             
include { MPOX as OTHER }         from './subworkflows/mpx'
include { SARSCOV2 }             from './subworkflows/sarscov2'

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
    cecret(
        ch_reads,
        ch_nanopore,
        ch_reference,
        ch_primer
    )
    ch_versions = ch_versions.mix(cecret.out.versions)
    ch_consensus = cecret.out.consensus.mix(ch_fastas).mix(ch_multifasta)

    qc(ch_reads,
      cecret.out.clean_reads,
      ch_kraken2_db,
      cecret.out.sam,
      cecret.out.trim_bam,
      ch_reference,
      ch_gff,
      ch_amplicon,
      ch_primer)

    ch_for_multiqc = cecret.out.for_multiqc.mix(qc.out.for_multiqc)
    ch_for_summary = qc.out.for_summary
    ch_versions    = ch_versions.mix(qc.out.versions)

    if ( params.species == 'sarscov2' ) {
      sarscov2(
        ch_consensus,
        cecret.out.trim_bam,
        ch_reference,
        ch_nextclade_dataset,
        ch_scripts
        )

      ch_for_multiqc = ch_for_multiqc.mix(sarscov2.out.for_multiqc)
      ch_for_summary = ch_for_summary.mix(sarscov2.out.for_summary)
      ch_versions    = ch_versions.mix(sarscov2.out.versions)
    
    } else if ( params.species == 'mpx') {
      mpx(
        ch_consensus, 
        ch_nextclade_dataset
        )
      
      ch_for_multiqc = ch_for_multiqc.mix(mpx.out.for_multiqc)
      ch_for_summary = ch_for_summary.mix(mpx.out.for_summary)
      ch_versions    = ch_versions.mix(mpx.out.versions)

    } else if ( params.species == 'other') {
      other(ch_consensus, 
        ch_nextclade_dataset
        )
      
      ch_for_multiqc = ch_for_multiqc.mix(other.out.for_multiqc)
      ch_for_summary = ch_for_summary.mix(other.out.for_summary)
      ch_versions    = ch_versions.mix(other.out.versions)

    } 

    if ( params.relatedness ) { 
      
      msa(
        ch_consensus,
        ch_reference_genome
        ) 

      tree      = msa.out.tree
      alignment = msa.out.msa
      matrix    = msa.out.matrix

      ch_for_multiqc = ch_for_multiqc.mix(msa.out.for_multiqc)
      ch_versions    = ch_versions.mix(msa.out.versions)

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
        ch_version_script
    )
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)

    summary(
      ch_for_summary
        .mix(fasta_prep.out.fastas).mix(cecret.out.consensus).collect().map{it -> tuple([it])}
        .combine(ch_script)
        .combine(ch_versions).collect().map{it -> tuple([it])}
        .combine(MULTIQC.out.files.ifEmpty([]).map{it -> tuple([it])})
    )

  emit:
    bam       = cecret.out.trim_bam
    consensus = ch_consensus
    tree      = tree
    alignment = alignment
    matrix    = matrix
}