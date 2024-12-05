include { FREYJA                        } from '../../modules/local/freyja'  
include { FREYJA_AGGREGATE as AGGREGATE } from '../../modules/local/freyja'  
include { PANGOLIN                      } from '../../modules/local/pangolin'  
include { PANGO_ALIASOR                 } from '../../modules/local/pango_aliasor' 
include { NEXTCLADE                     } from '../../modules/local/nextclade' 
include { NEXTCLADE_DATASET as DATASET  } from '../../modules/local/nextclade' 
include { UNZIP                         } from '../../modules/local/local' 
include { VADR                          } from '../../modules/local/vadr'

workflow SARSCOV2 {
    take:
        ch_fastas
        ch_bam
        ch_reference_genome
        ch_input_dataset
        ch_script

    main:
        ch_versions = Channel.empty()

        VADR(ch_fastas.collect())
        ch_versions = ch_versions.mix(VADR.out.versions)

        PANGOLIN(ch_fastas.collect())
        ch_versions = ch_versions.mix(PANGOLIN.out.versions)

        PANGO_ALIASOR(pangolin.out.pangolin_file)
        ch_versions = ch_versions.mix(PANGO_ALIASOR.out.versions)
        
        if ( params.download_nextclade_dataset ) {
            DATASET()
            ch_dataset = DATASET.out.dataset
            ch_versions = ch_versions.mix(DATASET.out.versions)
        } else {
            UNZIP(ch_input_dataset)
            ch_dataset = UNZIP.out.dataset
        }
        
        NEXTCLADE(ch_fastas.collect(), ch_dataset)
        ch_versions = ch_versions.mix(NEXTCLADE.out.versions)

        FREYJA(ch_bam.map{it -> tuple(it[0], it[1])}.combine(ch_reference_genome))
        ch_versions = ch_versions.mix(FREYJA.out.versions.first())

        AGGREGATE(FREYJA.out.demix.collect(), ch_script)
        ch_versions = ch_versions.mix(AGGREGATE.out.versions)

    emit:
        dataset     = ch_dataset
        for_multiqc = PANGOLIN.out.pangolin_file.mix(NEXTCLADE.out.nextclade_file).mix(AGGREGATE.out.for_multiqc)
        for_summary = AGGREGATE.out.aggregated_freyja_file.mix(VADR.out.vadr_file).mix(PANGO_ALIASOR.out.results).mix(NEXTCLADE.out.nextclade_file).mix(PANGOLIN.out.pangolin_file)
        versions    = ch_versions
}
