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

include { paramsHelp; validateParameters; paramsSummaryLog } from 'plugin/nf-schema'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow {

    main:
    if (params.help) {
        log.info paramsHelp("nextflow run UPHL-BioNGS/Cecret -profile docker --sample_sheet samplesheet.csv --outdir cecret")
        exit 0
    }

    // Validate parameters and print the summary log
    validateParameters()
    log.info paramsSummaryLog(workflow)


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
    if (workflow.success) {
        log.info """
        =============================================================================
        Cecret pipeline execution completed successfully
        =============================================================================
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}

        Primary Output Locations:
        -----------------------------------------------------------------------------
        • Results Summary (CSV)  : ${params.outdir}/cecret_results.csv
        • MultiQC Report         : ${params.outdir}/multiqc/multiqc_report.html
        • Consensus Genomes      : ${params.outdir}/consensus/
        • Pathogen Typing        : ${params.outdir}/pangolin/ | ${params.outdir}/nextclade/
        • Wastewater Abundance   : ${params.outdir}/freyja/
        • Phylogeny              : ${params.outdir}/iqtree/
        -----------------------------------------------------------------------------
        """
    } else {
        log.info """
        =============================================================================
        Cecret pipeline execution failed
        =============================================================================
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        
        Error message:
        ${workflow.errorMessage ?: 'No specific error message provided by Nextflow.'}

        Check the '.nextflow.log' file in the execution directory for full details.
        -----------------------------------------------------------------------------
        """
    }
}