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
        ch_vadr(fastas.collect())
        ch_pangolin(fastas.collect())
        ch_nextclade(fastas.collect())
        ch_freyja(bam.combine(reference_genome))
        ch_freyja_aggregate(freyja.out.freyja_demix.collect())

    emit:
        pangolin_file   = pangolin.out.pangolin_file
        nextclade_file  = nextclade.out.nextclade_file
        vadr_file       = vadr.out.vadr_file
        dataset         = nextclade.out.prepped_nextalign
        freyja_file     = freyja_aggregate.out.aggregated_freyja_file
}



  