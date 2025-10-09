//
// Subworkflow with functionality specific to the Cecret pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREP } from '../../../modules/local/local'
include { TEST } from '../../../subworkflows/local/test'


def determine_type(it) { 
    if (it == 'single') { 
        return true 
    } else { 
        return false 
    } 
}


// Define all presets and their corresponding files
def primer_presets = [
    'midnight_idt_V1' : [
        reference: "${workflow.projectDir}/genomes/MN908947.3.fasta",
        gff: "${workflow.projectDir}/genomes/MN908947.3.gff",
        primer_bed :  "${workflow.projectDir}/schema/midnight_idt_V1_SARS-CoV-2.primer.bed",
        amplicon_bed: "${workflow.projectDir}/schema/midnight_idt_V1_SARS-CoV-2.insert.bed",
    ],
    'midnight_ont_V1'  : [
        reference: "${workflow.projectDir}/genomes/MN908947.3.fasta",
        gff: "${workflow.projectDir}/genomes/MN908947.3.gff",
        primer_bed :  "${workflow.projectDir}/schema/midnight_ont_V1_SARS-CoV-2.primer.bed",
        amplicon_bed: "${workflow.projectDir}/schema/midnight_ont_V1_SARS-CoV-2.insert.bed",
    ],
    'midnight_ont_V2'  : [
        reference: "${workflow.projectDir}/genomes/MN908947.3.fasta",
        gff: "${workflow.projectDir}/genomes/MN908947.3.gff",
        primer_bed :  "${workflow.projectDir}/schema/midnight_ont_V2_SARS-CoV-2.primer.bed",
        amplicon_bed: "${workflow.projectDir}/schema/midnight_ont_V2_SARS-CoV-2.insert.bed",
    ],
    'midnight_ont_V3'  : [
        reference: "${workflow.projectDir}/genomes/MN908947.3.fasta",
        gff: "${workflow.projectDir}/genomes/MN908947.3.gff",
        primer_bed :  "${workflow.projectDir}/schema/midnight_ont_V3_SARS-CoV-2.primer.bed",
        amplicon_bed: "${workflow.projectDir}/schema/midnight_ont_V3_SARS-CoV-2.insert.bed",
    ],
    'ncov_V3'  : [
        reference: "${workflow.projectDir}/genomes/MN908947.3.fasta",
        gff: "${workflow.projectDir}/genomes/MN908947.3.gff",
        primer_bed :  "${workflow.projectDir}/schema/ncov_V3_nCoV-2019.primer.bed",
        amplicon_bed: "${workflow.projectDir}/schema/ncov_V3_nCoV-2019.insert.bed",
    ],
    'ncov_V4'  : [
        reference: "${workflow.projectDir}/genomes/MN908947.3.fasta",
        gff: "${workflow.projectDir}/genomes/MN908947.3.gff",
        primer_bed :  "${workflow.projectDir}/schema/ncov_V4_SARS-CoV-2.primer.bed",
        amplicon_bed: "${workflow.projectDir}/schema/ncov_V4_SARS-CoV-2.insert.bed",
    ],
    'ncov_V4.1'  : [
        reference: "${workflow.projectDir}/genomes/MN908947.3.fasta",
        gff: "${workflow.projectDir}/genomes/MN908947.3.gff",
        primer_bed :  "${workflow.projectDir}/schema/ncov_V4.1_SARS-CoV-2.primer.bed",
        amplicon_bed: "${workflow.projectDir}/schema/ncov_V4.1_SARS-CoV-2.insert.bed",
    ],
    'ncov_V5.3.2' : [
        reference: "${workflow.projectDir}/genomes/MN908947.3.fasta",
        gff: "${workflow.projectDir}/genomes/MN908947.3.gff",
        primer_bed :  "${workflow.projectDir}/schema/ncov_V5.3.2_SARS-CoV-2.primer.bed",
        amplicon_bed: "${workflow.projectDir}/schema/ncov_V5.3.2_SARS-CoV-2.insert.bed",
    ],
    'mpx_primalseq' : [
        reference: "${workflow.projectDir}/genomes/NC_063383.1.fasta",
        gff: "${workflow.projectDir}/genomes/NC_063383.1.gff3",
        primer_bed :  "${workflow.projectDir}/schema/mpx_primalseq_primer.bed",
        amplicon_bed: "${workflow.projectDir}/schema/mpx_primalseq_insert.bed",
    ],
    'mpx_yale' : [
        reference: "${workflow.projectDir}/genomes/MT903345.1.fasta",
        gff: "${workflow.projectDir}/genomes/MT903345.1.gff",
        primer_bed :  "${workflow.projectDir}/schema/mpx_yale_primer.bed",
        amplicon_bed: "${workflow.projectDir}/schema/mpx_yale_insert.bed",
    ],
    'mpx_idt' : [
        reference: "${workflow.projectDir}/genomes/NC_063383.1.fasta",
        gff: "${workflow.projectDir}/genomes/NC_063383.1.gff3",
        primer_bed :  "${workflow.projectDir}/schema/mpx_idt_primer.bed",
        amplicon_bed: "${workflow.projectDir}/schema/mpx_idt_insert.bed",
    ]
]

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALIZE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow INITIALIZE {
    main:
    ch_versions = Channel.empty()

    //# For aesthetics - and, yes, we are aware that there are better ways to write 
    //# this than a bunch of 'println' statements
    println('  ____ _____ ____ ____  _____ _____')
    println(' / ___| ____/ ___|  _ \\| ____|_   _|')
    println('| |   |  _|| |   | |_) |  _|   | |')
    println('| |___| |__| |___|  _ <| |___  | |')
    println(' \\____|_____\\____|_| \\_\\_____| |_|')

    println("Version: ${workflow.manifest.version}")
    println('')
    println('Currently using the Cecret workflow for use with corresponding reference genome.\n')
    println('Author: Erin Young')
    println('email: eriny@utah.gov')
    println('')

    println('Cecret is named after a real lake!')
    println('Visit https://www.alltrails.com/trail/us/utah/cecret-lake-trail to learn more.')
    println('Not everyone can visit in person, so here is some ASCII art of a mountain lake.')
    println(' ')
    println(' ')
    println('             /\\')
    println('            /  \\      /\\        /\\')
    println('           /    \\    /  \\      /  \\  /\\')
    println('          /      \\  /    \\    /    \\/  \\')
    println('   /\\    /   /\\   \\/      \\  /      \\   \\')
    println('  /  \\  /   /  \\           \\/        \\   \\')
    println(' /    \\/   /    \\    ~~~~~~~~         \\   \\')
    println('/         /      \\   ~~~~~~~~~~~      /\\   \\')
    println('            ~~~~~~~~   ~~~~~~~~~~    /  \\   \\')
    println('   ~~~~~~~~~~~~~~~~     ~~~~~~~~    /    \\   \\')
    println('  ~~~~~~~~~~~~~~~~       ~~~~~~    /      \\   \\')
    println('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    println('    ~~~~~~~~      ~~~~~~~~~~~~     ~~~~~~~~~~')
    println('     ~~~~~~~~~~~~   ~~~~~~~~~~ ~~~~~~~~  ~~~~')
    println(' ')

    //
    // Warn about deprecated params
    //

    //# Warning people about legacy params for a few versions. This was put here 3.6.20230418
    params.kraken2_organism                     = false
    if (params.kraken2_organism ) {
        println('WARNING : params.kraken2_organism no longer does anything!')
    }
    params.bedtools_multicov                    = false
    if (params.bedtools_multicov ) {
        println('WARNING : params.bedtools_multicov no longer does anything!')
    }
    params.bedtools_multicov_options            = false
    if (params.bedtools_multicov_options ) {
        println('WARNING : params.bedtools_multicov_options no longer does anything!')
    }
    params.freyja_boot_options                  = false
    if (params.freyja_boot_options  ) {
        println('WARNING : params.freyja_boot_options no longer does anything!')
    }
    params.nextalign_options                    = false
    if (params.nextalign_options ) {
        println('WARNING : params.nextalign_options no longer does anything!')
    }
    params.pango_collapse_options               = false
    if (params.pango_collapse_options ) {
        println('WARNING : params.pango_collapse_options no longer does anything!')
        println('WARNING : pango_collapse has been replaced by pango_aliasor')
    }
    params.medcpus                              = false
    if (params.medcpus ) {
        println('WARNING : params.medcpus no longer does anything!')
        println('WARNING : adjust computational resources with a config file')
    }
    params.maxcpus                              = false 
    if (params.maxcpus ) { 
        println('WARNING : params.maxcpus no longer does anything!')
        println('WARNING : adjust computational resources with a config file') 
    }
    if (params.msa == "nextalign" ) {
        println('WARNING : setting params.msa to nextalign no longer does anything!')
        println('WARNING : Use params.msa == "nextclade" instead!')
    }

    //
    // Create channels for input files
    //

    //# getting input files
    if ( params.sample_sheet ) { 
        Channel
            .fromPath(
                params.sample_sheet, 
                type: 'file', 
                checkIfExists: true
            )
            .view { "Sample sheet found : ${it}" }
            .splitCsv( header: true, sep: ',' )
            .map { it -> 
                def meta = [id:it.sample, single_end:determine_type(it.fastq_2)]
                tuple(meta, [ file("${it.fastq_1}", checkIfExists: true), file("${it.fastq_2}") ]) }
            .branch {
                single     : it[1] =~ /single/
                multifasta : it[1] =~ /multifasta/
                fasta      : it[1] =~ /fasta/
                ont        : it[1] =~ /ont/       
                paired     : true 
            }
            .set { inputs }

        ch_paired_reads = inputs.paired.map{     it -> tuple(it[0], it[1])}
        ch_single_reads = inputs.single.map{     it -> tuple(it[0], it[1][0])}
        ch_fastas       = inputs.fasta.map{      it -> tuple(it[0], it[1][0])}
        ch_nanopore     = inputs.ont.map{        it -> tuple(it[0], it[1][0])}
        ch_multifastas  = inputs.multifasta.map{ it -> tuple(it[1][0])}

    } else {

        //# Ensuring that reads and single_reads are not set to the same directory
        if ( params.reads && params.reads == params.single_reads ) {
            println("'params.reads' and 'params.single_reads' cannot point to the same directory!")
            println("'params.reads' is set to ${params.reads}")
            println("'params.single_reads' is set to ${params.single_reads}")
            exit 1
        }

        if ( params.fastas && params.fastas == params.multifastas ) {
            println("'params.fastas' and 'params.multifastas' cannot point to the same directory!")
            println("'params.fastas' is set to ${params.fastas}")
            println("'params.multifastas' is set to ${params.multifastas}")
            exit 1
        }


        // looking for paired-end fastq files
        if (params.reads) {
            Channel
                .fromFilePairs(["${params.reads}/*_R{1,2}*.{fastq,fastq.gz,fq,fq.gz}",
                        "${params.reads}/*_{1,2}*.{fastq,fastq.gz,fq,fq.gz}"], size: 2 )
                .unique()
                .map { it ->
                    def meta = [id:it[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), single_end:false] 
                    tuple( meta, [
                        file(it[1][0], checkIfExists: true), 
                        file(it[1][1], checkIfExists: true)]
                    )
                }
                .unique()
                .ifEmpty{
                    println('FATAL : No input files were found!')
                    println("No paired-end fastq files were found at ${params.reads}. Set 'params.reads' to directory with paired-end reads")
                    exit 1
                }
                .view { "Paired-end fastq files found : ${it[0].id}" }
                .set { ch_paired_reads }
        } else {
            ch_paired_reads = Channel.empty()
        }

        // looking for paired-end fastq files
        if (params.single_reads) {
            Channel
                .fromPath("${params.single_reads}/*.{fastq,fastq.gz,fq,fq.gz}")
                .map { it -> 
                    def meta = [id:it.simpleName, single_end:true] 
                    tuple( meta, file(it, checkIfExists: true))
                }
                .unique()
                .ifEmpty{
                    println('FATAL : No input files were found!')
                    println("No single-end fastq files were found at ${params.single_reads}. Set 'params.single_reads' to directory with single-end reads")
                    exit 1
                }
                .view { "Single-end fastq files found : ${it[0].id}" }
                .set { ch_single_reads }
        } else {
            ch_single_reads = Channel.empty()
        }

        // looking for nanopore fastq files
        if (params.nanopore) {
            Channel
                .fromPath("${params.nanopore}/*.{fastq,fastq.gz,fq,fq.gz}")
                .map { it -> 
                    def meta = [id:it.simpleName, single_end:true] 
                    tuple( meta, file(it, checkIfExists: true))
                }
                .unique()
                .ifEmpty{
                    println('FATAL : No input files were found!')
                    println("No nanopore fastq files were found at ${params.nanopore}. Set 'params.nanopore' to directory with nanopore reads")
                    exit 1
                }
                .view { "Nanopore fastq files found : ${it[0].id}" }
                .set { ch_nanopore }

        } else {
            ch_nanopore = Channel.empty()
        }

        // looking for fasta files
        if (params.fastas) {
            Channel
                .fromPath("${params.fastas}/*{.fa,.fasta,.fna}", type:'file')
                .map { it -> 
                    def meta = [id:it.simpleName, single_end:null] 
                    tuple( meta, file(it, checkIfExists: true))
                }
                .unique()
                .ifEmpty{
                    println('FATAL : No input files were found!')
                    println("No fasta files were found at ${params.fastas}. Set 'params.fastas' to directory with fastas.")
                    exit 1
                }
                .view { "Fasta files found : ${it[0].id}" }
                .set { ch_fastas }
        } else {
            ch_fastas = Channel.empty()
        }

        // looking for multifasta files
        if (params.multifastas) {
            Channel
                .fromPath("${params.multifastas}/*{.fa,.fasta,.fna}", type:'file')
                .unique()
                .ifEmpty{
                    println('FATAL : No input files were found!')
                    println("No multifasta files were found at ${params.multifastas}. Set 'params.multifastas' to directory with multifastas.")
                    exit 1
                }
                .view { "Multifasta files found : ${it}" }
                .set { ch_multifastas }
        } else {
            ch_multifastas = Channel.empty()
        }
    }

    // getting fastq files from sra accessions (Illumina only)
    ch_sra_accessions = Channel.from( params.sra_accessions )

    // getting genome accessions (ncbi virus is default)
    ch_genome_accessions = Channel.from( params.genome_accessions )

    //# getting a reference genome file
    if (params.reference_genome){
        Channel
            .fromPath(
                params.reference_genome, 
                type:'file'
            )
            .ifEmpty{
                println("No reference genome was selected. Set with 'params.reference_genome'")
                exit 1
            }
            .set { ch_reference }
    } else {
        if ( params.primer_set && primer_presets.containsKey(params.primer_set)) {
            Channel
                .fromPath(
                    primer_presets[params.primer_set].reference, 
                    type:'file', 
                    checkIfExists: true
                )
                .set { ch_reference }
        } else {
            println("WARN: No reference genome was selected. Set with 'params.reference_genome'")
            println("Or set species to one with an included genome ('sarscov2' or 'mpx')")
            ch_reference = Channel.empty()
        } 
    }
    ch_reference.view { "Reference Genome : $it"}

    //# getting the gff file for ivar variants
    if ( params.ivar_variants ) {
        if (params.gff) {
            Channel
            .fromPath(
                params.gff, 
                type:'file'
            )
            .ifEmpty{
                println("No gff file was selected. Set with 'params.reference_genome'")
                exit 1
            }
            .set { ch_gff }

        } else {
            if ( params.primer_set && primer_presets.containsKey(params.primer_set)) {
                Channel
                    .fromPath(
                        primer_presets[params.primer_set].gff, 
                        type:'file', 
                        checkIfExists: true
                    )
                    .set { ch_gff }
            } else {
                println("No gff file was selected. Set with 'params.gff'")
                println("Or set 'params.species' to one with an included genome ('sarscov2' or 'mpx')")
                println("Or bypass this message completely by setting 'params.ivar_variants = false'")
                exit 1

                // so nextflow doesn't throw an error too
                ch_gff = Channel.empty()
            } 
        }
        ch_gff.view { "GFF file : $it"}
    } else {
        ch_gff = Channel.empty()
    }

    if ( params.trimmer != 'none' ) {
        //# Getting the primer file
        if (params.primer_bed) {
            Channel.fromPath(
                params.primer_bed, 
                type:'file'
            )
            .ifEmpty{
                println("A bedfile for primers is required. Set with 'params.primer_bed'.")
                println("or use a provided primer set in ${primer_presets.keySet()}")
                exit 1
            }
            .set { ch_primer_bed }
        } else if ( params.primer_set && primer_presets.containsKey(params.primer_set)) {
            Channel
                .fromPath(
                    primer_presets[params.primer_set].primer_bed, 
                    type:'file', 
                    checkIfExists: true
                )
                .set { ch_primer_bed }
        } else {
            println("WARN: No primers were found!")
            println("Set primer schema with 'params.primer_bed' or specify to 'none' if primers were not used")
            println("Or use included primer set by setting 'params.primer_set' to one of ${primer_presets.keySet()}")
        
            ch_primer_bed = Channel.empty()
        }
        ch_primer_bed.view { "Primer BedFile : $it"}
    } else {
        ch_primer_bed = Channel.empty()
    }

    //# Getting the amplicon bedfile
    if ( params.aci ) {
        if (params.amplicon_bed) {
            Channel
                .fromPath(
                    params.amplicon_bed, 
                    type:'file'
                )
                .ifEmpty{
                    println("A bedfile for amplicons is required. Set with 'params.amplicon_bed'.")
                    println("Or set params.aci = false to skip this.")
                    exit 1
                }
                .set { ch_amplicon_bed } 
        } else if ( params.primer_set && primer_presets.containsKey(params.primer_set)) {
            Channel
                .fromPath(
                    primer_presets[params.primer_set].amplicon_bed, 
                    type:'file', 
                    checkIfExists: true
                )
                .set { ch_amplicon_bed } 
        } else {
            println("An amplicon bedfile wasn't found!")
            println("Set amplicon schema with 'params.amplicon_bed'")
            println("Or use included primer set by setting 'params.primer_set' to one of ${primer_presets.keySet()}")
            println("Or set params.aci = false to skip this.")
            exit 1

            // so nextflow doesn't throw an error too
            ch_amplicon_bed = Channel.empty()
        }
        ch_amplicon_bed.view { "Amplicon BedFile : $it"}
    } else {
        ch_amplicon_bed = Channel.empty()
    }

    //# scripts for legacy reasons
    scripts = [
        "${workflow.projectDir}/bin/combine_results.py",
        "${workflow.projectDir}/bin/freyja_graphs.py",
        "${workflow.projectDir}/bin/versions.py"
    ]

    ch_scripts = Channel.fromPath(scripts, type: 'file').collect()

    if ( params.kraken2_db ) {
        Channel
            .fromPath(params.kraken2_db, type:'dir')
            .view { "Kraken2 database : ${it}" }
            .set{ ch_kraken2_db }
    } else {
        ch_kraken2_db = Channel.empty()
    }

    if ( params.predownloaded_nextclade_dataset ) {
        Channel
            .fromPath(
                params.predownloaded_nextclade_dataset, 
                type:'file'
            )
            .ifEmpty{
                println("Dataset file could not be found at ${params.predownloaded_nextclade_dataset}.")
                println("Please set nextclade dataset file with 'params.predownloaded_nextclade_dataset'")
                exit 1
            }
            .view { "NextClade Dataset File : ${it}" }
            .set { ch_nextclade_dataset }
    } else if (! params.download_nextclade_dataset ) {
        Channel
            .fromPath(
                "${workflow.projectDir}/data/sars.zip", 
                type:'file',
                checkIfExists: true
            )
            .view { "NextClade Dataset File : ${it}" }
            .set { ch_nextclade_dataset }
    } else {
        ch_nextclade_dataset = Channel.empty()
    }

    ch_paired_reads
        .mix(ch_single_reads)
        .unique()
        .set { ch_reads }

    if (params.sample_sheet || params.fastas ) {
        PREP(ch_fastas)
        ch_versions = ch_versions.mix(PREP.out.versions)
        ch_prepped_fastas = PREP.out.fastas
    } else {
        ch_prepped_fastas = Channel.empty()
    }

    if ( ! params.sra_accessions.isEmpty() || ! params.genome_accessions.isEmpty() ) { 
        TEST(
            ch_sra_accessions,
            ch_genome_accessions
        )

        ch_reads = ch_reads.mix(TEST.out.reads)
        ch_prepped_fastas = ch_prepped_fastas.mix(TEST.out.fastas)
        ch_versions = ch_versions.mix(TEST.out.versions)
    }

    emit:
    reads               = ch_reads // channel: [meta, reads]
    fastas              = ch_prepped_fastas // channel: [meta, fasta]
    multifasta          = ch_multifastas // channel: fasta
    nanopore            = ch_nanopore // channel: [meta, reads]
    reference           = ch_reference // channel: fasta
    gff                 = ch_gff // channel: gff file
    primer              = ch_primer_bed // channel: bedfile
    amplicon            = ch_amplicon_bed // channel: bedfile 
    versions            = ch_versions // channel: value
    kraken2_db          = ch_kraken2_db // channel: path
    scripts             = ch_scripts // channel: [scripts]
    nextclade_dataset   = ch_nextclade_dataset // channel: file
}
