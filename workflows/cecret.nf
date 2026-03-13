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
include { OTHER         } from '../subworkflows/local/other'                             
include { OTHER as MPOX } from '../subworkflows/local/other'
include { SARSCOV2      } from '../subworkflows/local/sarscov2'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CECRET {

    take:
    ch_reads // channel: [meta, reads]
    ch_fastas // channel: fasta
    ch_multifasta // channel: fasta
    ch_nanopore // channel: [meta, reads]
    ch_reference // channel: fasta
    ch_gff // channel: gff file
    ch_primer // channel: bedfile
    ch_amplicon // channel: bedfile
    ch_versions // channel: value
    ch_kraken2_db // channel: path
    ch_scripts // channel: [scripts]
    ch_nextclade_dataset // channel: file

    main:
    log.info """

Running the complete Cecret Pipeline. This master workflow orchestrates consensus 
generation, quality control, species-specific typing, and phylogenetic analysis.

More information can be found at https://github.com/UPHL-BioNGS/Cecret/wiki/

If you run into issues or have feature requests, please contact the developers or 
submit an issue on GitHub at https://github.com/UPHL-BioNGS/Cecret/issues

┏━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ subworkflow        ┃ description                                                       ┃
┣━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
┃ CONSENSUS          ┃ Generates consensus sequences from raw reads (Illumina/Nanopore). ┃
┃ QC                 ┃ Performs comprehensive quality control and variant calling.       ┃
┃ SARSCOV2 / MPOX    ┃ Runs species-specific typing, validation, and abundance estimates.┃
┃ OTHER              ┃ Handles clade and lineage assignment for alternative pathogens.   ┃
┃ MSA                ┃ Aligns consensus sequences and builds a phylogenetic tree.        ┃
┃ MULTIQC            ┃ Aggregates QC metrics and logs into a single HTML report.         ┃
┃ SUMMARY            ┃ Compiles all pipeline outputs into a final results spreadsheet.   ┃
┗━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

"""
    ch_versions    = ch_versions
    ch_for_multiqc = channel.empty()
    ch_for_summary = channel.empty()
    ch_bam         = channel.empty()
    ch_consensus   = ch_fastas.mix(ch_multifasta)
    ch_wo_mltifna  = ch_fastas

    if (params.reads || params.single_reads || params.sample_sheet || params.nanopore || ! params.sra_accessions.isEmpty()) {
      CONSENSUS(
          ch_reads,
          ch_nanopore,
          ch_reference,
          ch_primer
      )
      ch_versions    = ch_versions.mix(CONSENSUS.out.versions)
      ch_consensus   = ch_consensus.mix(CONSENSUS.out.consensus)
      ch_for_multiqc = ch_for_multiqc.mix(CONSENSUS.out.for_multiqc)
      ch_wo_mltifna  = ch_wo_mltifna.mix(CONSENSUS.out.consensus)

      QC(ch_reads,
        CONSENSUS.out.clean_reads,
        ch_kraken2_db.ifEmpty([]),
        CONSENSUS.out.just_bam,
        CONSENSUS.out.trim_bam,
        ch_reference,
        ch_gff,
        ch_amplicon,
        ch_primer)

      ch_for_multiqc = ch_for_multiqc.mix(QC.out.for_multiqc)
      ch_for_summary = ch_for_summary.mix(QC.out.for_summary)
      ch_versions    = ch_versions.mix(QC.out.versions)

      ch_bam         = CONSENSUS.out.trim_bam
    } 

    if ( params.species == 'sarscov2' ) {
      SARSCOV2(
        ch_consensus,
        ch_bam,
        ch_reference,
        ch_nextclade_dataset,
        ch_scripts
        )

      ch_for_multiqc = ch_for_multiqc.mix(SARSCOV2.out.for_multiqc)
      ch_for_summary = ch_for_summary.mix(SARSCOV2.out.for_summary)
      ch_versions    = ch_versions.mix(SARSCOV2.out.versions)
    
    } else if ( params.species == 'mpx') {
      MPOX(
        ch_consensus,
        ch_bam,
        ch_reference,
        ch_nextclade_dataset,
        ch_scripts
        )
      
      ch_for_multiqc = ch_for_multiqc.mix(MPOX.out.for_multiqc)
      ch_for_summary = ch_for_summary.mix(MPOX.out.for_summary)
      ch_versions    = ch_versions.mix(MPOX.out.versions)

    } else if ( (params.species != 'sarscov2' && params.species != 'mpx' ) && (params.nextclade_dataset != 'sarscov2' || params.vadr_reference != 'sarscov2')) {
      OTHER(
        ch_consensus,
        ch_bam,
        ch_reference,
        ch_nextclade_dataset,
        ch_scripts
        )
      
      ch_for_multiqc = ch_for_multiqc.mix(OTHER.out.for_multiqc)
      ch_for_summary = ch_for_summary.mix(OTHER.out.for_summary)
      ch_versions    = ch_versions.mix(OTHER.out.versions)

    } 

    if ( params.relatedness ) { 

      MSA(
        ch_consensus,
        ch_reference.ifEmpty([])
        ) 

      tree      = MSA.out.tree
      alignment = MSA.out.msa
      matrix    = MSA.out.matrix

      ch_for_multiqc = ch_for_multiqc.mix(MSA.out.for_multiqc)
      ch_versions    = ch_versions.mix(MSA.out.versions)

    } else {
      tree      = channel.empty()
      alignment = channel.empty()
      matrix    = channel.empty()
    }

    ch_versions
      .collectFile(
        keepHeader: false,
        name: "collated_versions.yml")
      .set { ch_collated_versions }

    if (params.multiqc){
      MULTIQC(
          ch_for_multiqc.mix(ch_collated_versions).collect(), 
          ch_scripts
      )
      ch_versions    = ch_versions.mix(MULTIQC.out.versions)

      SUMMARY(
        ch_for_summary.mix(ch_wo_mltifna).collect().map{it -> tuple([it])}
          .combine(ch_scripts.map{it -> tuple([it])})
          .combine(MULTIQC.out.files.ifEmpty([]).map{it -> tuple([it])})
      )
    }

  emit:
    bam       = ch_bam
    consensus = ch_consensus
    tree      = tree
    alignment = alignment
    matrix    = matrix
}
