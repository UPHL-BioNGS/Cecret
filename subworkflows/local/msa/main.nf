include { MAFFT }       from '../../../modules/local/mafft'
include { HEATCLUSTER } from '../../../modules/local/heatcluster'
include { IQTREE }      from '../../../modules/local/iqtree' 
include { PHYTREEVIZ }  from '../../../modules/local/phytreeviz'
include { SNPDISTS }    from '../../../modules/local/snp-dists'

workflow MSA {
  take:
    ch_fasta // channel: fasta
    ch_reference_genome // channel: fasta
    
  main:
    log.info """

Running relatedness and phylogeny analysis. This workflow takes the generated 
consensus sequences, aligns them to the reference, and calculates genetic 
distances and phylogenetic trees.

Relevant params and their values:
- 'params.iqtree_outgroup' : ${params.iqtree_outgroup}
    - Explicitly designate the outgroup for IQTREE.

┏━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ process            ┃ description                                                       ┃
┣━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
┃ MAFFT              ┃ Aligns consensus sequences to create a multiple sequence alignment┃
┃ IQTREE             ┃ Infers a maximum-likelihood phylogenetic tree from the alignment. ┃
┃ PHYTREEVIZ         ┃ Generates a visual rendering of the generated phylogenetic tree.  ┃
┃ SNP-DISTS          ┃ Calculates a pairwise SNP distance matrix from the alignment.     ┃
┃ HEATCLUSTER        ┃ Creates a clustered heatmap from the SNP distance matrix.         ┃
┗━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

"""

    ch_versions = channel.empty()
    ch_nwk      = channel.empty()
    ch_matrix   = channel.empty()
    ch_multiqc  = channel.empty()

    // use mafft for msa
    // there used to be nextalign for an option
    // TODO : add in USHER as an option
    // TODO : allow the user to use nextclades aligned fasta

    ch_fasta
      .collect()
      .filter { it -> it.size() >= 3 }
      .set{ch_collected_fastas }

    if ( params.msa == 'mafft' ) {
      MAFFT(ch_collected_fastas, ch_reference_genome)
      ch_msa = MAFFT.out.msa
      ch_versions = ch_versions.mix(MAFFT.out.versions)
    } else {
      ch_msa = channel.empty()
    }
    
    // run iqtree
    if (params.iqtree) {
      IQTREE(ch_msa)
      ch_versions = ch_versions.mix(IQTREE.out.versions)
      ch_nwk      = ch_nwk.mix(IQTREE.out.newick)
    }

    if (params.phytreeviz) {
      // run phytreeviz for a basic visualization
      PHYTREEVIZ(ch_nwk)
      ch_versions = ch_versions.mix(PHYTREEVIZ.out.versions)
      ch_multiqc  = ch_multiqc.mix(PHYTREEVIZ.out.for_multiqc)
    }

    if (params.snpdists) {
      // run snp-dists for a snp matrix
      SNPDISTS(ch_msa)
      ch_matrix = ch_matrix.mix(SNPDISTS.out.matrix)
      ch_versions = ch_versions.mix(SNPDISTS.out.versions)
    }

    if (params.heatcluster && params.snpdists ) {
      // run heatcluster for a basic visualization
      HEATCLUSTER(SNPDISTS.out.matrix)
      ch_versions = ch_versions.mix(HEATCLUSTER.out.versions)
      ch_multiqc  = ch_multiqc.mix(HEATCLUSTER.out.for_multiqc)
    }

  emit:
    tree        = ch_nwk
    matrix      = ch_matrix
    msa         = ch_msa
    for_multiqc = ch_multiqc
    versions    = ch_versions
}
