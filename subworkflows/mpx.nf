include { nextclade }                    from '../modules/nextclade' addParams(params)
include { nextclade_dataset as dataset } from '../modules/nextclade' addParams(params)
include { unzip }                        from '../modules/cecret'    addParams(params)
include { vadr }                         from '../modules/vadr'      addParams(params)

workflow mpx {
  take:
    ch_fastas
    ch_nextclade_dataset

  main:
    vadr(ch_fastas.collect())
    ch_versions = vadr.out.versions

    if ( params.download_nextclade_dataset ) {
      dataset()
      ch_dataset = dataset.out.dataset
      ch_versions = ch_versions.mix(dataset.out.versions)
    } else {
      unzip(ch_input_dataset)
      ch_dataset = unzip.out.dataset
    }
        
    nextclade(ch_fastas.collect(), ch_dataset)
    ch_versions = ch_versions.mix(nextclade.out.versions)

  emit:
    for_multiqc = nextclade.out.nextclade_file
    for_summary = vadr.out.vadr_file
    prealigned  = nextclade.out.prealigned
    dataset     = ch_dataset
    versions    = ch_versions
}