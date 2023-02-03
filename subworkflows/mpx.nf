include { vadr }      from '../modules/vadr'      addParams(params)
include { nextclade } from '../modules/nextclade' addParams(params)

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