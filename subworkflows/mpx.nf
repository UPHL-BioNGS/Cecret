include { vadr }      from '../modules/vadr'      addParams(params)
include { nextclade } from '../modules/nextclade' addParams(params)

workflow mpx {
  take:
    ch_fastas

  main:
    vadr(ch_fastas.collect())
    nextclade(ch_fastas.collect())

  emit:
    for_multiqc = nextclade.out.nextclade_file
    for_summary = nextclade.out.nextclade_file.mix(vadr.out.vadr_file)
    dataset     = nextclade.out.prepped_nextalign
}