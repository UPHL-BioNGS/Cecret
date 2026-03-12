include { FREYJA                        } from '../../../modules/local/freyja'  
include { FREYJA_AGGREGATE as AGGREGATE } from '../../../modules/local/freyja'  
include { PANGOLIN                      } from '../../../modules/local/pangolin'  
include { PANGO_ALIASOR                 } from '../../../modules/local/pango_aliasor' 
include { NEXTCLADE                     } from '../../../modules/local/nextclade' 
include { NEXTCLADE_DATASET as DATASET  } from '../../../modules/local/nextclade' 
include { UNZIP                         } from '../../../modules/local/local' 
include { VADR                          } from '../../../modules/local/vadr'

workflow SARSCOV2 {
    take:
        ch_fastas // channel: fasta
        ch_bam // channel: [meta, bam]
        ch_reference_genome // channel: fasta
        ch_input_dataset // channel: zipped file
        ch_script // channel: workflow scripts

    main:

    log.info """

Relevant params and their values:
- 'params.minimum_reads' : ${params.minimum_reads}
    - Any samples with fewer than this will not be included in other steps.


Running SARS-CoV-2 specific analysis. This workflow performs lineage assignment, 
clade typing, sequence validation, and wastewater abundance estimation tailored 
for SARS-CoV-2 genomes.

┏━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ process            ┃ description                                                       ┃
┣━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
┃ VADR               ┃ Validates sequences against a SARS-CoV-2 reference model.         ┃
┃ PANGOLIN           ┃ Assigns PANGO lineages to the generated consensus sequences.      ┃
┃ PANGO_ALIASOR      ┃ Translates PANGO aliases (e.g., BA.1) to full lineage strings.    ┃
┃ DATASET / UNZIP    ┃ Fetches or extracts the Nextclade SARS-CoV-2 dataset.             ┃
┃ NEXTCLADE          ┃ Assigns Nextstrain clades, calls mutations, and performs QC.      ┃
┃ FREYJA             ┃ Estimates relative lineage abundance from mixed BAM files.        ┃
┃ AGGREGATE          ┃ Aggregates Freyja abundance outputs across multiple samples.      ┃
┗━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

"""
        ch_versions = channel.empty()
        ch_multiqc  = channel.empty()
        ch_summary  = channel.empty()

        if (params.vadr) {
            // running vadr
            VADR(ch_fastas.collect())
            ch_versions = ch_versions.mix(VADR.out.versions)
            ch_summary  = ch_summary.mix(VADR.out.vadr_file)
        }

        if (params.pangolin) {
            // running pangolin
            PANGOLIN(ch_fastas.collect())
            ch_versions = ch_versions.mix(PANGOLIN.out.versions)
            ch_summary  = ch_summary.mix(PANGOLIN.out.pangolin_file)
            ch_multiqc  = ch_multiqc.mix(PANGOLIN.out.pangolin_file)
        }

        if (params.pango_aliasor && params.pangolin) {
            // running pango aliasor
            PANGO_ALIASOR(PANGOLIN.out.pangolin_file)
            ch_versions = ch_versions.mix(PANGO_ALIASOR.out.versions)
            ch_summary  = ch_summary.mix(PANGO_ALIASOR.out.results)
        }

        // running nextclade
        if (params.nextclade) {
            if ( params.download_nextclade_dataset && !params.predownloaded_nextclade_dataset ) {
                DATASET()
                ch_dataset = DATASET.out.dataset
                ch_versions = ch_versions.mix(DATASET.out.versions)
            } else {
                UNZIP(ch_input_dataset)
                ch_dataset = UNZIP.out.dataset
                ch_versions = ch_versions.mix(UNZIP.out.versions)
            }
            
            NEXTCLADE(ch_fastas.collect(), ch_dataset)
            ch_versions = ch_versions.mix(NEXTCLADE.out.versions)
            ch_summary  = ch_summary.mix(NEXTCLADE.out.nextclade_file)
            ch_multiqc  = ch_multiqc.mix(NEXTCLADE.out.nextclade_file)
        }

        if (params.freyja) {
            // running freyja
            FREYJA(ch_bam.map{it -> tuple(it[0], it[1])}.combine(ch_reference_genome))
            ch_versions = ch_versions.mix(FREYJA.out.versions.first())

            if (params.freyja_aggregate) {
                AGGREGATE(FREYJA.out.demix.collect(), ch_script)
                ch_versions = ch_versions.mix(AGGREGATE.out.versions)
                ch_multiqc  = ch_multiqc.mix(AGGREGATE.out.for_multiqc)
                ch_summary  = ch_summary.mix(AGGREGATE.out.aggregated_freyja_file)
            }
        }

    emit:
        for_multiqc = ch_multiqc
        for_summary = ch_summary
        versions    = ch_versions
}
