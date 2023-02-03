include { mafft }     from '../modules/mafft'     addParams(params)
include { nextalign } from '../modules/nextalign' addParams(params)
include { iqtree2 }   from '../modules/iqtree2'   addParams(params)
include { snpdists }  from '../modules/snp-dists' addParams(params)

workflow msa {
  take:
  fasta
  reference_genome
  dataset
  
  main:
  if ( params.msa == 'nextalign' ) {
    nextalign(fasta.collect(), dataset)
    msa = nextalign.out.msa
  } else if ( params.msa == 'mafft' ) {
    mafft(fasta.collect(), reference_genome)
    msa = mafft.out.msa
  }
  iqtree2(msa)
  snpdists(msa)

  emit:
  tree   = iqtree2.out.tree
  matrix = snpdists.out.matrix
  msa    = msa
}
