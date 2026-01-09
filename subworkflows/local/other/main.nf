include { FREYJA_AGGREGATE as AGGREGATE } from '../../../modules/local/freyja'
include { FREYJA_PATHOGEN as FREYJA     } from '../../../modules/local/freyja'
include { FREYJA_UPDATE                 } from '../../../modules/local/freyja'
include { NEXTCLADE                     } from '../../../modules/local/nextclade'
include { NEXTCLADE_DATASET as DATASET  } from '../../../modules/local/nextclade'
include { UNZIP                         } from '../../../modules/local/local' 
include { VADR                          } from '../../../modules/local/vadr'

workflow OTHER {
  take:
  ch_fastas // channel: fasta
  ch_bam // channel: [meta, bam]
  ch_reference_genome // channel: fasta
  ch_input_dataset // channel: zipped file
  ch_script // channel: workflow scripts

  main:
  // create some empty channels for optional results
  ch_versions    = Channel.empty()
  ch_for_summary = Channel.empty()
  ch_for_multiqc = Channel.empty()

  // run vadr only if vadr is expected to run
  if (params.vadr_reference && params.vadr_reference != 'sarscov2') {
    VADR(ch_fastas.collect())
    ch_versions    = VADR.out.versions
    ch_for_summary = ch_for_summary.mix(VADR.out.vadr_file)
  }

  // run freyja only if freyja is expected to run
  if (params.freyja_pathogen && params.freyja_pathogen != 'SARS-CoV-2') {
    // running freyja
    // TODO : FIX WHEN NOT BROKEN
    // FREYJA_UPDATE()
    // ch_versions = ch_versions.mix(FREYJA_UPDATE.out.versions)

    //FREYJA(ch_bam.map{it -> tuple(it[0], it[1])}.combine(ch_reference_genome).combine(FREYJA_UPDATE.out.db))
    FREYJA(ch_bam.map{it -> tuple(it[0], it[1])}.combine(ch_reference_genome))
    ch_versions = ch_versions.mix(FREYJA.out.versions.first())

    if (params.freyja_aggregate) {
      AGGREGATE(FREYJA.out.demix.collect(), ch_script)
      ch_versions    = ch_versions.mix(AGGREGATE.out.versions)
      ch_for_multiqc = ch_for_multiqc.mix(AGGREGATE.out.for_multiqc)
      ch_for_summary = ch_for_summary.mix(AGGREGATE.out.aggregated_freyja_file)
    }
  }

  // run nextclade only if nextclade is expected to run
  if (( params.nextclade_dataset && params.nextclade_dataset != 'sars-cov-2' && params.download_nextclade_dataset ) || params.predownloaded_nextclade_dataset ) {

    // running nextclade
    if ( params.download_nextclade_dataset && !params.predownloaded_nextclade_dataset ) {
      DATASET()
      ch_dataset  = DATASET.out.dataset
      ch_versions = ch_versions.mix(DATASET.out.versions)
    } else {
      UNZIP(ch_input_dataset)
      ch_dataset  = UNZIP.out.dataset
      ch_versions = ch_versions.mix(UNZIP.out.versions)
    }
    
    NEXTCLADE(ch_fastas.collect(), DATASET.out.dataset)
    ch_versions    = ch_versions.mix(NEXTCLADE.out.versions)
    ch_for_multiqc = ch_for_multiqc.mix(NEXTCLADE.out.nextclade_file)
    ch_for_summary = ch_for_summary.mix(NEXTCLADE.out.nextclade_file)
  }

  emit:
    for_multiqc = ch_for_multiqc
    for_summary = ch_for_summary
    versions    = ch_versions
}