include { NEXTCLADE }                    from '../../modules/local/nextclade'
include { NEXTCLADE_DATASET as DATASET } from '../../modules/local/nextclade'
include { VADR }                         from '../../modules/local/vadr'

workflow OTHER {
  take:
  ch_fastas // channel: fasta 
  // ch_input_dataset // future goal of zipped directory of nextclade dataset

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

  // run nextclade only if nextclade is expected to run
  // TODO : Allow user to download and zip their own nextclade dataset
  if ( params.nextclade_dataset && params.nextclade_dataset != 'sars-cov-2' && params.download_nextclade_dataset ) {

    DATASET()
    ch_versions = ch_versions.mix(DATASET.out.versions)
    
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