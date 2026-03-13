include { FREYJA                        } from '../../../modules/local/freyja'
include { FREYJA_AGGREGATE              } from '../../../modules/local/freyja'
include { FREYJA_PATHOGEN               } from '../../../modules/local/freyja'
include { FREYJA_UPDATE                 } from '../../../modules/local/freyja'
include { NEXTCLADE                     } from '../../../modules/local/nextclade'
include { NEXTCLADE_DATASET as DATASET  } from '../../../modules/local/nextclade'
include { UNZIP                         } from '../../../modules/local/local' 
include { VADR                          } from '../../../modules/local/vadr'

workflow OTHER {
  take:
  ch_fastas // channel: fasta
  ch_bam // channel: [meta, bam]
  ch_reference_genome // channel: fasta
  ch_input_dataset // channel: zipped file
  ch_script // channel: workflow scripts

  main:
  log.info """

Running specific species or custom pathogen analysis. This workflow performs
sequence validation, clade assignment, and lineage abundance estimation for 
organisms other than the default SARS-CoV-2.

Relevant params and their values:

- 'params.species' : ${params.species}
    - Designates subworkflows
- 'params.nextclade_dataset': ${params.nextclade_dataset}
    - Designate which dataset to download in NEXTCLADE_DATASET.
    - will be ignored if value is set to 'sars-cov-2'.
    - See Nextclade documentation at 
      https://docs.nextstrain.org/projects/nextclade/en/stable/user/datasets.html to see 
      available datasets.
- 'params.freyja_pathogen': ${params.freyja_pathogen}
    - Designate which pathogen to download in FREYJA_UPDATE.
    - Will be ignored if value is set to 'SARS-CoV-2'.
    - See Freyja's documentation at 
      https://andersen-lab.github.io/Freyja/src/usage/update.html to see available 
      pathogens.
- 'params.vadr_reference': ${params.vadr_reference}
    - Designate with reference to use in the VADR process and corresponding container.
    - Will be ignored if value is equal to 'sarscov2'.
    - See available images at https://hub.docker.com/r/staphb/vadr/tags.
- 'params.download_nextclade_dataset' : ${params.download_nextclade_dataset}
    - Will used nextclade to download the dataset according to 'params.nextclade_dataset' 
      in NEXTCLADE_DATASET.
- 'params.predownloaded_nextclade_dataset' : ${params.predownloaded_nextclade_dataset}
    - Allows the user to use an existing nextclade dataset in a zipped directory.
    - See https://github.com/UPHL-BioNGS/Cecret/wiki/Usage#nextclade-datasets for more 
      information.

┏━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ process            ┃ description                                                       ┃
┣━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
┃ VADR               ┃ Validates viral sequences and annotates expected features/errors. ┃
┃ NEXTCLADE_DATASET  ┃ Downloads the requested Nextclade dataset for the target pathogen.┃
┃ NEXTCLADE          ┃ Performs clade assignment, mutation calling, and QC.              ┃
┃ FREYJA_UPDATE      ┃ Updates the Freyja database for the specified pathogen.           ┃
┃ FREYJA             ┃ Estimates lineage abundances from BAM files (e.g., wastewater).   ┃
┃ FREYJA_AGGREGATE   ┃ Aggregates Freyja abundance outputs across multiple samples.      ┃
┗━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

"""

  // create some empty channels for optional results
  ch_versions    = channel.empty()
  ch_for_summary = channel.empty()
  ch_for_multiqc = channel.empty()

  // run vadr only if vadr is expected to run
  if (params.vadr_reference && params.vadr_reference != 'sarscov2') {
    VADR(ch_fastas.collect())
    ch_versions    = VADR.out.versions
    ch_for_summary = ch_for_summary.mix(VADR.out.vadr_file)
  }

  // run freyja only if freyja is expected to run
  if (params.freyja_pathogen && params.freyja_pathogen != 'SARS-CoV-2') {

    if (params.freyja_update) {
      FREYJA_UPDATE(params.freyja_pathogen)
      ch_versions = ch_versions.mix(FREYJA_UPDATE.out.versions)

      FREYJA_PATHOGEN(ch_bam.map{it -> tuple(it[0], it[1])}.combine(ch_reference_genome).combine(FREYJA_UPDATE.out.db))
      ch_versions = ch_versions.mix(FREYJA_PATHOGEN.out.versions.first())
      ch_freyja_out = FREYJA_PATHOGEN.out.demix
    } else {
      FREYJA(ch_bam.map{it -> tuple(it[0], it[1])}.combine(ch_reference_genome))
      ch_versions = ch_versions.mix(FREYJA.out.versions.first())
      ch_freyja_out = FREYJA.out.demix
    }

    if (params.freyja_aggregate) {
      FREYJA_AGGREGATE(ch_freyja_out.collect(), ch_script)
      ch_versions    = ch_versions.mix(FREYJA_AGGREGATE.out.versions)
      ch_for_multiqc = ch_for_multiqc.mix(FREYJA_AGGREGATE.out.for_multiqc)
      ch_for_summary = ch_for_summary.mix(FREYJA_AGGREGATE.out.aggregated_freyja_file)
    }
  }

  // run nextclade only if nextclade is expected to run
  if (( params.nextclade_dataset && params.nextclade_dataset != 'sars-cov-2' && params.download_nextclade_dataset ) || params.predownloaded_nextclade_dataset ) {

    // running nextclade
    if ( params.download_nextclade_dataset && !params.predownloaded_nextclade_dataset ) {
      DATASET()
      ch_dataset  = DATASET.out.dataset
      ch_versions = ch_versions.mix(DATASET.out.versions)
    } else {
      UNZIP(ch_input_dataset)
      ch_dataset  = UNZIP.out.dataset
      ch_versions = ch_versions.mix(UNZIP.out.versions)
    }
    
    NEXTCLADE(ch_fastas.collect(), ch_dataset)
    ch_versions    = ch_versions.mix(NEXTCLADE.out.versions)
    ch_for_multiqc = ch_for_multiqc.mix(NEXTCLADE.out.nextclade_file)
    ch_for_summary = ch_for_summary.mix(NEXTCLADE.out.nextclade_file)
  }

  emit:
    for_multiqc = ch_for_multiqc
    for_summary = ch_for_summary
    versions    = ch_versions
}
