include { vadr }      from '../modules/vadr'      addParams(vadr: params.vadr, vadr_options: params.vadr_options, vadr_reference: params.vadr_reference, vadr_mdir: params.vadr_mdir)
include { pangolin }  from '../modules/pangolin'  addParams(pangolin: params.pangolin, pangolin_options: params.pangolin_options)
include { nextclade } from '../modules/nextclade' addParams(nextclade: params.nextclade, nextclade_options: params.nextclade_options, nextclade_dataset: params.nextclade_dataset)

workflow annotation {
  take:
  fastas

  main:
  vadr(fastas.collect())
  pangolin(fastas.collect())
  nextclade(fastas.collect())

  emit:
  pangolin_file   = pangolin.out.pangolin_file
  nextclade_file  = nextclade.out.nextclade_file
  vadr_file       = vadr.out.vadr_file
  dataset         = nextclade.out.prepped_nextalign
}
