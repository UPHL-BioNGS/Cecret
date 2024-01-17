include { mafft }       from '../modules/mafft'       addParams(params)
include { heatcluster } from '../modules/heatcluster' addParams(params)
include { iqtree2 }     from '../modules/iqtree2'     addParams(params)
include { phytreeviz }  from '../modules/phytreeviz'  addParams(params)
include { snpdists }    from '../modules/snp-dists'   addParams(params)

workflow msa {
  take:
    ch_fasta
    ch_reference_genome
    ch_prealigned
    
  main:
    ch_msa = Channel.empty()
    ch_nwk = Channel.empty()

    if ( params.msa == 'nextclade' ) {
      ch_msa = ch_prealigned.map { it -> tuple(it[0])}
      ch_nwk = ch_prealigned.map { it -> tuple(it[1])}
    } else if ( params.msa == 'mafft' ) {
      mafft(ch_fasta.collect(), ch_reference_genome)
      ch_msa = mafft.out.msa
      iqtree2(ch_msa)
      ch_nwk = iqtree2.out.newick
    }
    
    phytreeviz(ch_nwk)
    snpdists(ch_msa)
    heatcluster(snpdists.out.matrix)

  emit:
    tree        = ch_nwk
    matrix      = snpdists.out.matrix
    msa         = ch_msa
    for_multiqc = phytreeviz.out.for_multiqc.mix(heatcluster.out.for_multiqc)
}
