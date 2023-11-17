include { mafft }       from '../modules/mafft'       addParams(params)
include { nextalign }   from '../modules/nextalign'   addParams(params)
include { heatcluster } from '../modules/heatcluster' addParams(params)
include { iqtree2 }     from '../modules/iqtree2'     addParams(params)
include { phytreeviz }  from '../modules/phytreeviz'  addParams(params)
include { snpdists }    from '../modules/snp-dists'   addParams(params)

workflow msa {
  take:
    ch_fasta
    ch_reference_genome
    ch_dataset
    
  main:
    if ( params.msa == 'nextalign' ) {
      nextalign(ch_fasta.collect(), ch_dataset)
      ch_msa = nextalign.out.msa
    } else if ( params.msa == 'mafft' ) {
      mafft(ch_fasta.collect(), ch_reference_genome)
      ch_msa = mafft.out.msa
    }
    iqtree2(ch_msa)
    phytreeviz(iqtree2.out.newick)
    snpdists(ch_msa)
    heatcluster(snpdists.out.matrix)

  emit:
    tree        = iqtree2.out.tree
    matrix      = snpdists.out.matrix
    msa         = ch_msa
    for_multiqc = phytreeviz.out.for_multiqc.mix(heatcluster.out.for_multiqc)
}
