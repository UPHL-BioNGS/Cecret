include { freyja_aggregate }             from '../modules/freyja'    addParams(params)
include { freyja_demix }                 from '../modules/freyja'    addParams(params)
include { freyja_variants }              from '../modules/freyja'    addParams(params)
include { pangolin }                     from '../modules/pangolin'  addParams(params)
include { nextclade }                    from '../modules/nextclade' addParams(params)
include { nextclade_dataset as dataset } from '../modules/nextclade' addParams(params)
include { unzip }                        from '../modules/cecret'    addParams(params)
include { vadr }                         from '../modules/vadr'      addParams(params)

workflow sarscov2 {
    take:
        ch_fastas
        ch_bam
        ch_reference_genome
        ch_input_dataset
        ch_freyja_script

    main:
        vadr(ch_fastas.collect())
        pangolin(ch_fastas.collect())
        
        if ( params.download_nextclade_dataset ) {
            dataset()
            ch_dataset = dataset.out.dataset
        } else {
            unzip(ch_input_dataset)
            ch_dataset = unzip.out.dataset
        }
        
        nextclade(ch_fastas.collect(), ch_dataset)
        freyja_variants(ch_bam.map{it -> tuple(it[0], it[1])}.combine(ch_reference_genome))
        freyja_demix(freyja_variants.out.variants)
        freyja_aggregate(freyja_demix.out.demix.collect(), ch_freyja_script)

    emit:
        dataset     = ch_dataset
        for_multiqc = pangolin.out.pangolin_file.mix(nextclade.out.nextclade_file).mix(freyja_aggregate.out.for_multiqc)
        for_summary = freyja_aggregate.out.aggregated_freyja_file.mix(vadr.out.vadr_file)
}
