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








/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALIZE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow INITIALIZE {
    main:
    ch_versions = channel.empty()

    log.info """
    =============================================================================
     ____ _____ ____ ____  _____ _____ 
    / ___| ____/ ___|  _ \\| ____|_   _|
   | |   |  _|| |   | |_) |  _|   | |  
   | |___| |__| |___|  _ <| |___  | |  
    \\____|_____\\____|_| \\_\\_____| |_|  
    =============================================================================
    Version: ${workflow.manifest.version}
    
    Currently using the Cecret workflow for use with corresponding reference genome.
    
    Author: Erin Young
    email: eriny@utah.gov
    
    Cecret is named after a real lake!
    Visit https://www.alltrails.com/trail/us/utah/cecret-lake-trail to learn more.
    Not everyone can visit in person, so here is some ASCII art of a mountain lake.
     
     
                 /\\
                /  \\      /\\        /\\
               /    \\    /  \\      /  \\  /\\
              /      \\  /    \\    /    \\/  \\
       /\\    /   /\\   \\/      \\  /      \\   \\
      /  \\  /   /  \\           \\/        \\   \\
     /    \\/   /    \\    ~~~~~~~~         \\   \\
    /         /      \\   ~~~~~~~~~~~      /\\   \\
                ~~~~~~~~   ~~~~~~~~~~    /  \\   \\
       ~~~~~~~~~~~~~~~~     ~~~~~~~~    /    \\   \\
      ~~~~~~~~~~~~~~~~       ~~~~~~    /      \\   \\
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ~~~~~~~~      ~~~~~~~~~~~~     ~~~~~~~~~~
         ~~~~~~~~~~~~   ~~~~~~~~~~ ~~~~~~~~  ~~~~
     
    """

    //
    // Warn about deprecated params
    //

    //# Warning people about legacy params for a few versions. This was put here 3.6.20230418
    params.kraken2_organism                     = false
    if (params.kraken2_organism ) {
        log.warn 'WARNING : params.kraken2_organism no longer does anything!'
    }
    params.bedtools_multicov                    = false
    if (params.bedtools_multicov ) {
        log.warn 'WARNING : params.bedtools_multicov no longer does anything!'
    }
    params.bedtools_multicov_options            = false
    if (params.bedtools_multicov_options ) {
        log.warn 'WARNING : params.bedtools_multicov_options no longer does anything!'
    }
    params.freyja_boot_options                  = false
    if (params.freyja_boot_options  ) {
        log.warn 'WARNING : params.freyja_boot_options no longer does anything!'
    }
    params.nextalign_options                    = false
    if (params.nextalign_options ) {
        log.warn 'WARNING : params.nextalign_options no longer does anything!'
    }
    params.pango_collapse_options               = false
    if (params.pango_collapse_options ) {
        log.warn 'WARNING : params.pango_collapse_options no longer does anything!'
        log.warn 'WARNING : pango_collapse has been replaced by pango_aliasor'
    }
    params.medcpus                              = false
    if (params.medcpus ) {
        log.warn 'WARNING : params.medcpus no longer does anything!'
        log.warn 'WARNING : adjust computational resources with a config file'
    }
    params.maxcpus                              = false 
    if (params.maxcpus ) { 
        log.warn 'WARNING : params.maxcpus no longer does anything!'
        log.warn 'WARNING : adjust computational resources with a config file'
    }
    if (params.msa == "nextalign" ) {
        log.warn 'WARNING : setting params.msa to nextalign no longer does anything!'
    }

    //
    // Create channels for input files
    //

    log.info """
------------------------------------------------------
Initializing Sample Input Files
------------------------------------------------------
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
┃ param                             ┃ value                            
┣━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
┃ outdir                            ┃ ${params.outdir ?: 'None'}
┃ sample_sheet                      ┃ ${params.sample_sheet ?: 'None'}
┃ reads                             ┃ ${params.reads ?: 'None'}
┃ single_reads                      ┃ ${params.single_reads ?: 'None'}
┃ fastas                            ┃ ${params.fastas ?: 'None'}
┃ multifastas                       ┃ ${params.multifastas ?: 'None'}
┃ nanopore                          ┃ ${params.nanopore ?: 'None'}
┃ sra_accessions                    ┃ ${params.sra_accessions ?: 'None'}
┃ genome_accessions                 ┃ ${params.genome_accessions ?: 'None'}
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    """


    if ( params.sample_sheet ) { 
        channel
            .fromPath(
                params.sample_sheet, 
                type: 'file', 
                checkIfExists: true
            )
            .view { item -> "Sample sheet found : ${item}" }
            .splitCsv( header: true, sep: ',' )
            .map { it -> 
                def is_single = (it.fastq_2 == '' || it.fastq_2 == 'single') ? true : false
                def meta = [id:it.sample, single_end:is_single]
                def files = is_single ? [file("${it.fastq_1}", checkIfExists: true)] : [file("${it.fastq_1}", checkIfExists: true), file("${it.fastq_2}", checkIfExists: true)]
                
                tuple(meta, files) 
            }
            .branch { row ->
                single     : row[1] =~ /single/
                multifasta : row[1] =~ /multifasta/
                fasta      : row[1] =~ /fasta/
                ont        : row[1] =~ /ont/       
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
            log.warn "'params.reads' and 'params.single_reads' cannot point to the same directory!"
            log.warn "'params.reads' is set to ${params.reads}"
            log.warn "'params.single_reads' is set to ${params.single_reads}"
            exit 1
        }

        if ( params.fastas && params.fastas == params.multifastas ) {
            log.warn "'params.fastas' and 'params.multifastas' cannot point to the same directory!"
            log.warn "'params.fastas' is set to ${params.fastas}"
            log.warn "'params.multifastas' is set to ${params.multifastas}"
            exit 1
        }


        // looking for paired-end fastq files
        if (params.reads) {
            channel
                .fromFilePairs(["${params.reads}/*_R{1,2}*.{fastq,fastq.gz,fq,fq.gz}",
                        "${params.reads}/*_{1,2}*.{fastq,fastq.gz,fq,fq.gz}"], size: 2 )
                .map { it ->
                    def meta = [id:it[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), single_end:false] 
                    tuple( meta, [
                        file(it[1][0], checkIfExists: true), 
                        file(it[1][1], checkIfExists: true)]
                    )
                }
                .unique()
                .ifEmpty{
                    log.warn 'FATAL : No input files were found!'
                    log.warn "No paired-end fastq files were found at ${params.reads}. Set 'params.reads' to directory with paired-end reads"
                    exit 1
                }
                .view { item -> "Paired-end fastq files found : ${item[0].id}" }
                .set { ch_paired_reads }
        } else {
            ch_paired_reads = channel.empty()
        }

        // looking for paired-end fastq files
        if (params.single_reads) {
            channel
                .fromPath("${params.single_reads}/*.{fastq,fastq.gz,fq,fq.gz}")
                .map { it -> 
                    def meta = [id:it.simpleName, single_end:true] 
                    tuple( meta, file(it, checkIfExists: true))
                }
                .unique()
                .ifEmpty{
                    log.warn 'FATAL : No input files were found!'
                    log.warn "No single-end fastq files were found at ${params.single_reads}. Set 'params.single_reads' to directory with single-end reads"
                    exit 1
                }
                .view { it -> "Single-end fastq files found : ${it[0].id}" }
                .set { ch_single_reads }
        } else {
            ch_single_reads = channel.empty()
        }

        // looking for nanopore fastq files
        if (params.nanopore) {
            channel
                .fromPath("${params.nanopore}/*.{fastq,fastq.gz,fq,fq.gz}")
                .map { it -> 
                    def meta = [id:it.simpleName, single_end:true] 
                    tuple( meta, file(it, checkIfExists: true))
                }
                .unique()
                .ifEmpty{
                    log.warn 'FATAL : No input files were found!'
                    log.warn "No nanopore fastq files were found at ${params.nanopore}. Set 'params.nanopore' to directory with nanopore reads"
                    exit 1
                }
                .view { it -> "Nanopore fastq files found : ${it[0].id}" }
                .set { ch_nanopore }

        } else {
            ch_nanopore = channel.empty()
        }

        // looking for fasta files
        if (params.fastas) {
            channel
                .fromPath("${params.fastas}/*{.fa,.fasta,.fna}", type:'file')
                .map { it -> 
                    def meta = [id:it.simpleName, single_end:null] 
                    tuple( meta, file(it, checkIfExists: true))
                }
                .unique()
                .ifEmpty{
                    log.warn 'FATAL : No input files were found!'
                    log.warn "No fasta files were found at ${params.fastas}. Set 'params.fastas' to directory with fastas."
                    exit 1
                }
                .view { it -> "Fasta files found : ${it[0].id}" }
                .set { ch_fastas }
        } else {
            ch_fastas = channel.empty()
        }

        // looking for multifasta files
        if (params.multifastas) {
            channel
                .fromPath("${params.multifastas}/*{.fa,.fasta,.fna}", type:'file')
                .unique()
                .ifEmpty{
                    log.warn 'FATAL : No input files were found!'
                    log.warn "No multifasta files were found at ${params.multifastas}. Set 'params.multifastas' to directory with multifastas."
                    exit 1
                }
                .view { it -> "Multifasta files found : ${it}" }
                .set { ch_multifastas }
        } else {
            ch_multifastas = channel.empty()
        }
    }

    // getting fastq files from sra accessions (Illumina only)
    ch_sra_accessions = params.sra_accessions
        ? channel.from(params.sra_accessions)
        : channel.empty()

    // getting genome accessions (ncbi virus is default)
    ch_genome_accessions = params.genome_accessions
        ? channel.from(params.genome_accessions)
        : channel.empty()

    log.info """
------------------------------------------------------
Initializing References and Databases
------------------------------------------------------
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
┃ param                             ┃ value                            
┣━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
┃ species                           ┃ ${params.species}
┃ primer_set                        ┃ ${params.primer_set ?: 'None'}
┃ reference_genome                  ┃ ${params.reference_genome ?: 'Set with Primer Set'}
┃ gff                               ┃ ${params.gff ?: 'Set with Primer Set'}
┃ primer_bed                        ┃ ${params.primer_bed ?: 'Set with Primer Set'}
┃ amplicon_bed                      ┃ ${params.amplicon_bed ?: 'Set with Primer Set'}
┃ kraken2_db                        ┃ ${params.kraken2_db ?: 'Included in Container'}
┃ predownloaded_nextclade_dataset   ┃ ${params.predownloaded_nextclade_dataset ?: 'Included SARS-CoV-2'}
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""
    // Define all presets and their corresponding files
    def primer_presets = [
        'midnight_idt_V1' : [
            reference: "${projectDir}/genomes/MN908947.3.fasta",
            gff: "${projectDir}/genomes/MN908947.3.gff",
            primer_bed :  "${projectDir}/schema/midnight_idt_V1_SARS-CoV-2.primer.bed",
            amplicon_bed: "${projectDir}/schema/midnight_idt_V1_SARS-CoV-2.insert.bed",
        ],
        'midnight_ont_V1'  : [
            reference: "${projectDir}/genomes/MN908947.3.fasta",
            gff: "${projectDir}/genomes/MN908947.3.gff",
            primer_bed :  "${projectDir}/schema/midnight_ont_V1_SARS-CoV-2.primer.bed",
            amplicon_bed: "${projectDir}/schema/midnight_ont_V1_SARS-CoV-2.insert.bed",
        ],
        'midnight_ont_V2'  : [
            reference: "${projectDir}/genomes/MN908947.3.fasta",
            gff: "${projectDir}/genomes/MN908947.3.gff",
            primer_bed :  "${projectDir}/schema/midnight_ont_V2_SARS-CoV-2.primer.bed",
            amplicon_bed: "${projectDir}/schema/midnight_ont_V2_SARS-CoV-2.insert.bed",
        ],
        'midnight_ont_V3'  : [
            reference: "${projectDir}/genomes/MN908947.3.fasta",
            gff: "${projectDir}/genomes/MN908947.3.gff",
            primer_bed :  "${projectDir}/schema/midnight_ont_V3_SARS-CoV-2.primer.bed",
            amplicon_bed: "${projectDir}/schema/midnight_ont_V3_SARS-CoV-2.insert.bed",
        ],
        'ncov_V3'  : [
            reference: "${projectDir}/genomes/MN908947.3.fasta",
            gff: "${projectDir}/genomes/MN908947.3.gff",
            primer_bed :  "${projectDir}/schema/ncov_V3_nCoV-2019.primer.bed",
            amplicon_bed: "${projectDir}/schema/ncov_V3_nCoV-2019.insert.bed",
        ],
        'ncov_V4'  : [
            reference: "${projectDir}/genomes/MN908947.3.fasta",
            gff: "${projectDir}/genomes/MN908947.3.gff",
            primer_bed :  "${projectDir}/schema/ncov_V4_SARS-CoV-2.primer.bed",
            amplicon_bed: "${projectDir}/schema/ncov_V4_SARS-CoV-2.insert.bed",
        ],
        'ncov_V4.1'  : [
            reference: "${projectDir}/genomes/MN908947.3.fasta",
            gff: "${projectDir}/genomes/MN908947.3.gff",
            primer_bed :  "${projectDir}/schema/ncov_V4.1_SARS-CoV-2.primer.bed",
            amplicon_bed: "${projectDir}/schema/ncov_V4.1_SARS-CoV-2.insert.bed",
        ],
        'ncov_V5.3.2' : [
            reference: "${projectDir}/genomes/MN908947.3.fasta",
            gff: "${projectDir}/genomes/MN908947.3.gff",
            primer_bed :  "${projectDir}/schema/ncov_V5.3.2_SARS-CoV-2.primer.bed",
            amplicon_bed: "${projectDir}/schema/ncov_V5.3.2_SARS-CoV-2.insert.bed",
        ],
        'mpx_primalseq' : [
            reference: "${projectDir}/genomes/NC_063383.1.fasta",
            gff: "${projectDir}/genomes/NC_063383.1.gff3",
            primer_bed :  "${projectDir}/schema/mpx_primalseq_primer.bed",
            amplicon_bed: "${projectDir}/schema/mpx_primalseq_insert.bed",
        ],
        'mpx_yale' : [
            reference: "${projectDir}/genomes/MT903345.1.fasta",
            gff: "${projectDir}/genomes/MT903345.1.gff",
            primer_bed :  "${projectDir}/schema/mpx_yale_primer.bed",
            amplicon_bed: "${projectDir}/schema/mpx_yale_insert.bed",
        ],
        'mpx_idt' : [
            reference: "${projectDir}/genomes/NC_063383.1.fasta",
            gff: "${projectDir}/genomes/NC_063383.1.gff3",
            primer_bed :  "${projectDir}/schema/mpx_idt_primer.bed",
            amplicon_bed: "${projectDir}/schema/mpx_idt_insert.bed",
        ]
    ]
    //# getting a reference genome file
    if (params.reference_genome){
        channel
            .fromPath(
                params.reference_genome, 
                type:'file'
            )
            .ifEmpty{
                log.fatal "No reference genome was selected. Set with 'params.reference_genome'"
                exit 1
            }
            .set { ch_reference }
    } else {
        if ( params.primer_set && primer_presets.containsKey(params.primer_set)) {
            channel
                .fromPath(
                    primer_presets[params.primer_set].reference, 
                    type:'file', 
                    checkIfExists: true
                )
                .set { ch_reference }
        } else {
            log.warn "WARN: No reference genome was selected. Set with 'params.reference_genome'"
            log.warn "Or set species to one with an included genome ('sarscov2' or 'mpx')"
            ch_reference = channel.empty()
        } 
    }
    ch_reference.view { it -> "Reference Genome : $it"}

    //# getting the gff file for ivar variants
    if ( params.ivar_variants ) {
        if (params.gff) {
            channel
            .fromPath(
                params.gff, 
                type:'file'
            )
            .ifEmpty{
                log.warn "No gff file was selected. Set with 'params.gff'"
                exit 1
            }
            .set { ch_gff }

        } else {
            if ( params.primer_set && primer_presets.containsKey(params.primer_set)) {
                channel
                    .fromPath(
                        primer_presets[params.primer_set].gff, 
                        type:'file', 
                        checkIfExists: true
                    )
                    .set { ch_gff }
            } else {
                log.warn "No gff file was selected. Set with 'params.gff'"
                log.warn "Or set 'params.species' to one with an included genome ('sarscov2' or 'mpx')"
                log.warn "Or bypass this message completely by setting 'params.ivar_variants = false'"
                exit 1

                // so nextflow doesn't throw an error too
                ch_gff = channel.empty()
            } 
        }
        ch_gff.view { it -> "GFF file : $it"}
    } else {
        ch_gff = channel.empty()
    }

    if ( params.trimmer != 'none' ) {

        //# Getting the primer file
        if (params.primer_bed) {
            channel.fromPath(
                params.primer_bed,
                type:'file'
            )
            .ifEmpty{
                log.warn "A bedfile for primers is required. Set with 'params.primer_bed'."
                log.warn "or use a provided primer set in ${primer_presets.keySet()}"
                exit 1
            }
            .set { ch_primer_bed }

        } else if ( params.primer_set ) {

            if (primer_presets.containsKey(params.primer_set)) {
                channel
                    .fromPath(
                        primer_presets[params.primer_set].primer_bed,
                        type:'file',
                        checkIfExists: true
                    )
                    .set { ch_primer_bed }

            } else {
                log.warn "Primer set ${params.primer_set} not found!"
                log.warn "Set primer schema with 'params.primer_bed' or specify to 'none' if primers were not used"
                log.warn "Or use included primer set by setting 'params.primer_set' to one of ${primer_presets.keySet()}"
                exit 1
            }

        } else {
            ch_primer_bed = channel.empty()
        }

        ch_primer_bed.view { it -> "Primer BedFile : $it" }

    } else {
        ch_primer_bed = channel.empty()
    }

    //# Getting the amplicon bedfile
    if ( params.aci ) {
        if (params.amplicon_bed) {
            channel
                .fromPath(
                    params.amplicon_bed, 
                    type:'file',
                    checkIfExists: true
                )
                .ifEmpty{
                    log.warn "A bedfile for amplicons is required. Set with 'params.amplicon_bed'."
                    log.warn "Or set params.aci = false to skip this."
                    exit 1
                }
                .set { ch_amplicon_bed } 
        } else if ( params.primer_set && primer_presets.containsKey(params.primer_set)) {
            channel
                .fromPath(
                    primer_presets[params.primer_set].amplicon_bed, 
                    type:'file', 
                    checkIfExists: true
                )
                .set { ch_amplicon_bed } 
        } else {
            log.warn "An amplicon bedfile wasn't found!"
            log.warn "Set amplicon schema with 'params.amplicon_bed'"
            log.warn "Or use included primer set by setting 'params.primer_set' to one of ${primer_presets.keySet()}"
            log.warn "Or set params.aci = false to skip this."
            exit 1

            // so nextflow doesn't throw an error too
            ch_amplicon_bed = channel.empty()
        }
        ch_amplicon_bed.view { it -> "Amplicon BedFile : $it"}
    } else {
        ch_amplicon_bed = channel.empty()
    }

    //# scripts for legacy reasons
    scripts = [
        "${projectDir}/bin/combine_results.py",
        "${projectDir}/bin/freyja_graphs.py",
        "${projectDir}/bin/versions.py"
    ]

    ch_scripts = channel.fromPath(scripts, type: 'file').collect()

    if ( params.kraken2_db ) {
        channel
            .fromPath(params.kraken2_db, type:'dir')
            .view { it -> "Kraken2 database : ${it}" }
            .set{ ch_kraken2_db }
    } else {
        ch_kraken2_db = channel.empty()
    }

    if ( params.predownloaded_nextclade_dataset ) {
        channel
            .fromPath(
                params.predownloaded_nextclade_dataset, 
                type:'file'
            )
            .ifEmpty{
                log.warn "Dataset file could not be found at ${params.predownloaded_nextclade_dataset}."
                log.warn "Please set nextclade dataset file with 'params.predownloaded_nextclade_dataset'"
                exit 1
            }
            .view { it -> "NextClade Dataset File : ${it}" }
            .set { ch_nextclade_dataset }
    } else if (! params.download_nextclade_dataset ) {
        channel
            .fromPath(
                "${projectDir}/data/sars.zip", 
                type:'file',
                checkIfExists: true
            )
            .view { it -> "NextClade Dataset File : ${it}" }
            .set { ch_nextclade_dataset }
    } else {
        ch_nextclade_dataset = channel.empty()
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
        ch_prepped_fastas = channel.empty()
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

    log.info """
------------------------------------------------------
Initializing Species Values
------------------------------------------------------
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
┃ param                             ┃ value                            
┣━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
┃ species                           ┃ ${params.species ?: 'None'}
┃ vadr_reference                    ┃ ${params.vadr_reference ?: 'None'}
┃ nextclade_dataset                 ┃ ${params.nextclade_dataset ?: 'None'}
┃ freyja_pathogen                   ┃ ${params.freyja_pathogen ?: 'None'}
┃ iqtree_outgroup                   ┃ ${params.iqtree_outgroup ?: 'None'}
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    """

    log.info """
------------------------------------------------------
Initializing Subworkflow Toggles
------------------------------------------------------
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
┃ param                             ┃ value                            
┣━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
┃ trimmer                           ┃ ${params.trimmer ?: 'None'}
┃ cleaner                           ┃ ${params.cleaner ?: 'None'}
┃ aligner                           ┃ ${params.aligner ?: 'None'}
┃ msa                               ┃ ${params.msa ?: 'None'}
┃ relatedness                       ┃ ${params.relatedness ?: 'None'}
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    """

    if (params.relatedness) {
        log.info "'params.relatedness' is set to true. All input files will be put through the Multiple Sequence Alignment (MSA) subworkflow."
    } else {
        log.info "FYI: The MSA subworkflow is skipped by default. To compare isolates with reasonable top hits, set 'params.relatedness' to true."
    }

    //
    // Dynamically check for unused parameters based on toggles
    //
    def toggle_dependencies = [
        'kraken2':       ['kraken2_db', 'kraken2_options'],
        'ivar_variants': ['gff', 'ivar_variants_options'],
        'vadr':          ['vadr_reference', 'vadr_mdir', 'vadr_options'],
        'nextclade':     ['nextclade_dataset', 'download_nextclade_dataset', 'predownloaded_nextclade_dataset'],
        'freyja':        ['freyja_pathogen', 'freyja_aggregate', 'freyja_update'],
        'relatedness':   ['msa', 'iqtree', 'phytreeviz', 'snpdists', 'heatcluster'],
        'pangolin':      ['pango_aliasor', 'pangolin_options']
    ]
    
    def ignored_list = []
    
    toggle_dependencies.each { toggle, dependents ->
        if (!params[toggle]) {
            dependents.each { dep -> 
                ignored_list << "┃ ${dep.padRight(33)} ┃ Disabled by: ${toggle}"
            }
        }
    }

    if (ignored_list.size() > 0) {
        log.info """
------------------------------------------------------
Ignored Parameters (Due to Disabled Processes)
------------------------------------------------------
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
┃ ignored parameter                 ┃ reason             
┣━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
${ignored_list.join('\n')}
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        """
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
