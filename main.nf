#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    erinyoung/Cecret
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/erinyoung/Cecret
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INITIALIZE } from './subworkflows/local/initialize'
include { CECRET     } from './workflows/cecret'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow {

    main:
    //
    // SUBWORKFLOW: Initialize files and tasks
    //
    INITIALIZE ()

    //
    // WORKFLOW: Run main workflow
    //
    CECRET (
      INITIALIZE.out.reads, // channel: [id, reads]
      INITIALIZE.out.fastas, // channel: [id, fasta]
      INITIALIZE.out.multifasta, // channel: fasta
      INITIALIZE.out.nanopore, // channel: [id, reads]
      INITIALIZE.out.reference, // channel: fasta
      INITIALIZE.out.gff, // channel: gff file
      INITIALIZE.out.primer, // channel: bedfile
      INITIALIZE.out.amplicon, // channel: bedfile
      INITIALIZE.out.for_version, // channel: value
      INITIALIZE.out.versions, // channel: value
      INITIALIZE.out.kraken2_db, // channel: path
      INITIALIZE.out.scripts, // channel: [scripts]
      INITIALIZE.out.nextclade_dataset // channel: file
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    println("Pipeline completed at: $workflow.complete")
    println("A summary of results can be found in a comma-delimited file: ${params.outdir}/cecret_results.csv")
    println("A summary of results can be found in a tab-delimited file: ${params.outdir}/cecret_results.txt")
    println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
}