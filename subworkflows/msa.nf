include { mafft }     from '../modules/mafft'     addParams(mafft_options: params.mafft_options)
include { nextalign } from '../modules/nextalign' addParams(nextalign_options: params.nextalign_options)
include { iqtree2 }   from '../modules/iqtree2'   addParams(iqtree2: params.iqtree2, iqtree2_outgroup: params.iqtree2_outgroup)
include { snpdists }  from '../modules/snp-dists' addParams(snpdists: params.snpdists, snpdists_options: params.snpdists_options)

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
}
