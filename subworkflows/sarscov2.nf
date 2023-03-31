include { freyja; freyja_aggregate }    from '../modules/freyja'    addParams(params)
include { vadr }                        from '../modules/vadr'      addParams(params)
include { pangolin }                    from '../modules/pangolin'  addParams(params)
include { nextclade }                   from '../modules/nextclade' addParams(params)

workflow sarscov2 {
    take:
        ch_fastas
        ch_bam
        ch_reference_genome

    main:
        vadr(ch_fastas.collect())
        pangolin(ch_fastas.collect())
        nextclade(ch_fastas.collect())
        freyja(ch_bam.map{it -> tuple(it[0], it[1])}.combine(ch_reference_genome))
        freyja_aggregate(freyja.out.freyja_demix.collect())

    emit:
        dataset     = nextclade.out.prepped_nextalign
        for_multiqc = pangolin.out.pangolin_file.mix(nextclade.out.nextclade_file)
        for_summary = pangolin.out.pangolin_file.mix(nextclade.out.nextclade_file).mix(freyja_aggregate.out.aggregated_freyja_file).mix(vadr.out.vadr_file)
}



  