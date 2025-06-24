include { MAFFT }       from '../../../modules/local/mafft'
include { HEATCLUSTER } from '../../../modules/local/heatcluster'
include { IQTREE2 }     from '../../../modules/local/iqtree2' 
include { PHYTREEVIZ }  from '../../../modules/local/phytreeviz'
include { SNPDISTS }    from '../../../modules/local/snp-dists'

workflow MSA {
  take:
    ch_fasta // channel: fasta
    ch_reference_genome // channel: fasta
    
  main:
    ch_versions = Channel.empty()
    ch_nwk      = Channel.empty()

    // use mafft for msa
    // there used to be nextalign for an option
    // TODO : add in USHER as an option
    // TODO : allow the user to use nextclades aligned fasta
    if ( params.msa == 'mafft' ) {
      MAFFT(ch_fasta.collect(), ch_reference_genome)
      ch_msa = MAFFT.out.msa
      ch_versions = ch_versions.mix(MAFFT.out.versions)
    } else {
      ch_msa = Channel.empty()
    }
    
    // run iqtree2
    if (params.iqtree2) {
      IQTREE2(ch_msa)
      ch_versions = ch_versions.mix(IQTREE2.out.versions)
      ch_nwk      = ch_nwk.mix(IQTREE2.out.newick)
    }

    if (params.phytreeviz) {
      // run phytreeviz for a basic visualization
      PHYTREEVIZ(ch_nwk)
      ch_versions = ch_versions.mix(PHYTREEVIZ.out.versions)
    }

    if (params.snpdists) {
      // run snp-dists for a snp matrix
      SNPDISTS(ch_msa)
      ch_versions = ch_versions.mix(SNPDISTS.out.versions)
    }

    if (params.heatcluster && params.snpdists ) {
      // run heatcluster for a basic visualization
      HEATCLUSTER(SNPDISTS.out.matrix)
      ch_versions = ch_versions.mix(HEATCLUSTER.out.versions)
    }

  emit:
    tree        = IQTREE2.out.newick
    matrix      = SNPDISTS.out.matrix
    msa         = ch_msa
    for_multiqc = PHYTREEVIZ.out.for_multiqc.mix(HEATCLUSTER.out.for_multiqc)
    versions    = ch_versions
}
