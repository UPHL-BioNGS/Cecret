#!/usr/bin/env nextflow


/*
 * Generate pipeline help directly from nextflow_schema.json.
 */
def generateHelp(boolean showHidden = false, boolean fullHelp = false) {

    def schemaFile = new File("${projectDir}/nextflow_schema.json")

    if (!schemaFile.exists()) {
        throw new IllegalStateException(
            "Cannot find nextflow_schema.json at: ${schemaFile}"
        )
    }

    def schema = new groovy.json.JsonSlurper().parse(schemaFile)
    def output = []

    output << ""
    output << schema.title
    output << schema.description
    output << ""
    output << "USAGE:"
    output << "  nextflow run UPHL-BioNGS/Cecret -profile docker --sample_sheet samplesheet.csv --outdir cecret"
    output << ""
    output << "OPTIONS:"

    // Follow the schema's allOf references to retrieve parameter groups.
    schema.allOf.each { groupRef ->

        def groupName = groupRef['$ref'].tokenize('/').last()
        def group = schema['$defs'][groupName]

        if (!group) {
            return
        }

        output << ""
        output << group.title.toUpperCase()

        if (group.description) {
            output << group.description
        }

        def requiredParams = group.required ?: []

        group.properties.each { paramName, option ->

            // Skip hidden parameters unless explicitly requested.
            if (option.hidden && !showHidden) {
                return
            }

            def description = option.description ?: "No description available."

            output << ""
            output << String.format(
                "  --%-32s %s",
                paramName,
                description
            )

            def details = []

            if (requiredParams.contains(paramName)) {
                details << "Required"
            }

            if (option.type) {
                def type = option.type instanceof List ?
                    option.type.join(" or ") : option.type
                details << "Type: ${type}"
            }

            if (option.containsKey('default')) {
                details << "Default: ${option.default}"
            }

            if (option['enum']) {
                details << "Choices: ${option['enum'].join(', ')}"
            }

            if (details) {
                output << "    ${details.join(' | ')}"
            }

            // Include extended help text in full-help mode.
            if (fullHelp && option.help_text) {
                output << "    ${option.help_text}"
            }
        }
    }

    output << ""
    output << "Run with --help_full to display all parameters."
    output << ""

    return output.join('\n')
}


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
    if (params.help || params.help_full) {
        def helpText = generateHelp(
            params.show_hidden || params.help_full,
            params.help_full
        )

        log.info(helpText)
        exit 0
    }

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
    
    log.info """
        Cecret Output Locations:
        -----------------------------------------------------------------------------
        • Results Summary (CSV)  : ${params.outdir}/cecret_results.csv
        • MultiQC Report         : ${params.outdir}/multiqc/multiqc_report.html
        • Consensus Genomes      : ${params.outdir}/consensus/
        • Pathogen Typing        : ${params.outdir}/pangolin/ | ${params.outdir}/nextclade/
        • Wastewater Abundance   : ${params.outdir}/freyja/
        • Phylogeny              : ${params.outdir}/iqtree/
        -----------------------------------------------------------------------------
"""

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
