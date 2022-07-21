include { freyja; freyja_aggregate }    from '../modules/freyja'    addParams(freyja: params.freyja, freyja_variants_options: params.freyja_variants_options, freyja_demix_options: params.freyja_demix_options, freyja_aggregate: params.freyja_aggregate, freyja_aggregate_options: params.freyja_aggregate_options, freyja_plot_options: params.freyja_plot_options)
include { vadr }                        from '../modules/vadr'      addParams(vadr: params.vadr, vadr_options: params.vadr_options, vadr_reference: params.vadr_reference, vadr_mdir: params.vadr_mdir)
include { pangolin }                    from '../modules/pangolin'  addParams(pangolin: params.pangolin, pangolin_options: params.pangolin_options)
include { nextclade }                   from '../modules/nextclade' addParams(nextclade: params.nextclade, nextclade_options: params.nextclade_options, nextclade_dataset: params.nextclade_dataset)

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



  