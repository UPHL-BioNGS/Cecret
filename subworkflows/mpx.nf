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
            
    if ( params.download_nextclade_dataset ) {
      dataset()
      ch_dataset = dataset.out.dataset
    } else {
      unzip(ch_input_dataset)
      ch_dataset = unzip.out.dataset
    }
        
    nextclade(ch_fastas.collect(), ch_dataset)

  emit:
    for_multiqc = nextclade.out.nextclade_file
    for_summary = vadr.out.vadr_file
    dataset     = ch_dataset
}