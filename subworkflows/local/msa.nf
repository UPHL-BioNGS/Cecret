include { MAFFT }       from '../../modules/local/mafft'
include { HEATCLUSTER } from '../../modules/local/heatcluster'
include { IQTREE2 }     from '../../modules/local/iqtree2' 
include { PHYTREEVIZ }  from '../../modules/local/phytreeviz'
include { SNPDISTS }    from '../../modules/local/snp-dists'

workflow MSA {
  take:
    ch_fasta
    ch_reference_genome
    
  main:
    ch_versions = Channel.empty()

    if ( params.msa == 'mafft' ) {
      MAFFT(ch_fasta.collect(), ch_reference_genome)
      ch_msa = MAFFT.out.msa
      ch_versions = ch_versions.mix(MAFFT.out.versions)
    } else {
      ch_msa = Channel.empty()
    }
    
    IQTREE2(ch_msa)
    ch_versions = ch_versions.mix(IQTREE2.out.versions)

    PHYTREEVIZ(IQTREE2.out.newick)
    ch_versions = ch_versions.mix(PHYTREEVIZ.out.versions)

    SNPDISTS(ch_msa)
    ch_versions = ch_versions.mix(SNPDISTS.out.versions)

    HEATCLUSTER(SNPDISTS.out.matrix)
    ch_versions = ch_versions.mix(HEATCLUSTER.out.versions)

  emit:
    tree        = IQTREE2.out.newick
    matrix      = SNPDISTS.out.matrix
    msa         = ch_msa
    for_multiqc = PHYTREEVIZ.out.for_multiqc.mix(HEATCLUSTER.out.for_multiqc)
    versions    = ch_versions
}
