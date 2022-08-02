include { vadr }      from '../modules/vadr'      addParams(vadr: params.vadr, vadr_options: params.vadr_options, vadr_reference: params.vadr_reference, vadr_mdir: params.vadr_mdir)
include { nextclade } from '../modules/nextclade' addParams(nextclade: params.nextclade, nextclade_options: params.nextclade_options, nextclade_dataset: params.nextclade_dataset)

workflow mpx {
  take:
  fastas

  main:
  vadr(fastas.collect())
  nextclade(fastas.collect())

  emit:
  nextclade_file  = nextclade.out.nextclade_file
  vadr_file       = vadr.out.vadr_file
  dataset         = nextclade.out.prepped_nextalign
}