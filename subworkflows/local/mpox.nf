include { NEXTCLADE }                    from '../../modules/local/nextclade'
include { NEXTCLADE_DATASET as DATASET } from '../../modules/local/nextclade'
include { VADR }                         from '../../modules/local/vadr'

workflow MPOX {
  take:
    ch_fastas

  main:
    VADR(ch_fastas.collect())
    ch_versions = VADR.out.versions

    if ( params.download_nextclade_dataset ) {
      DATASET()
      ch_dataset = DATASET.out.dataset
      ch_versions = ch_versions.mix(DATASET.out.versions)
    } else {
      ch_dataset = Channel.empty()
    }
    
    NEXTCLADE(ch_fastas.collect(), ch_dataset)
    ch_versions = ch_versions.mix(NEXTCLADE.out.versions)

  emit:
    for_multiqc = NEXTCLADE.out.nextclade_file
    for_summary = VADR.out.vadr_file.mix(NEXTCLADE.out.nextclade_file)
    dataset     = ch_dataset
    versions    = ch_versions
}