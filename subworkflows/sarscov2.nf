include { freyja; freyja_aggregate }    from '../modules/freyja'    addParams(params)
include { vadr }                        from '../modules/vadr'      addParams(params)
include { pangolin }                    from '../modules/pangolin'  addParams(params)
include { nextclade }                   from '../modules/nextclade' addParams(params)

workflow sarscov2 {
    take:
    fastas
    bam
    reference_genome

    main:
    vadr(fastas.collect())
    pangolin(fastas.collect())
    nextclade(fastas.collect())
    freyja(bam.combine(reference_genome))
    freyja_aggregate(freyja.out.freyja_demix.collect())

    emit:
    pangolin_file   = pangolin.out.pangolin_file
    nextclade_file  = nextclade.out.nextclade_file
    vadr_file       = vadr.out.vadr_file
    dataset         = nextclade.out.prepped_nextalign
    freyja_file     = freyja_aggregate.out.aggregated_freyja_file
}



  