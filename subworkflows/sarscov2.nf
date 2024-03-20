include { freyja_aggregate }             from '../modules/freyja'        addParams(params)
include { freyja_demix }                 from '../modules/freyja'        addParams(params)
include { freyja_variants }              from '../modules/freyja'        addParams(params)
include { pangolin }                     from '../modules/pangolin'      addParams(params)
include { pango_collapse }               from '../modules/pangocollapse' addParams(params)
include { nextclade }                    from '../modules/nextclade'     addParams(params)
include { nextclade_dataset as dataset } from '../modules/nextclade'     addParams(params)
include { unzip }                        from '../modules/cecret'        addParams(params)
include { vadr }                         from '../modules/vadr'          addParams(params)

workflow sarscov2 {
    take:
        ch_fastas
        ch_bam
        ch_reference_genome
        ch_input_dataset
        ch_freyja_script

    main:
        ch_versions = Channel.empty()

        vadr(ch_fastas.collect())
        pangolin(ch_fastas.collect())
        pango_collapse(pangolin.out.pangolin_file)

        ch_versions = ch_versions.mix(vadr.out.versions)
        ch_versions = ch_versions.mix(pangolin.out.versions)
        ch_versions = ch_versions.mix(pango_collapse.out.versions)
        
        if ( params.download_nextclade_dataset ) {
            dataset()
            ch_dataset = dataset.out.dataset
            ch_versions = ch_versions.mix(dataset.out.versions)
        } else {
            unzip(ch_input_dataset)
            ch_dataset = unzip.out.dataset
        }
        
        nextclade(ch_fastas.collect(), ch_dataset)
        freyja_variants(ch_bam.map{it -> tuple(it[0], it[1])}.combine(ch_reference_genome))
        freyja_demix(freyja_variants.out.variants)
        freyja_aggregate(freyja_demix.out.demix.collect(), ch_freyja_script)
        
        ch_versions = ch_versions.mix(nextclade.out.versions)
        ch_versions = ch_versions.mix(freyja_variants.out.versions.first())
        ch_versions = ch_versions.mix(freyja_demix.out.versions.first())
        ch_versions = ch_versions.mix(freyja_aggregate.out.versions)

    emit:
        dataset     = ch_dataset
        prealigned  = nextclade.out.prealigned
        for_multiqc = pangolin.out.pangolin_file.mix(nextclade.out.nextclade_file).mix(freyja_aggregate.out.for_multiqc)
        for_summary = freyja_aggregate.out.aggregated_freyja_file.mix(vadr.out.vadr_file).mix(pango_collapse.out.results)
        versions    = ch_versions
}
