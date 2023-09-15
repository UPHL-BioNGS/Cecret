#!/usr/bin/env nextflow

//# For aesthetics - and, yes, we are aware that there are better ways to write this than a bunch of 'println' statements
println('') 
println('  ____ _____ ____ ____  _____ _____')
println(' / ___| ____/ ___|  _ \\| ____|_   _|')
println('| |   |  _|| |   | |_) |  _|   | |')
println('| |___| |__| |___|  _ <| |___  | |')
println(' \\____|_____\\____|_| \\_\\_____| |_|')

println('Version: ' + workflow.manifest.version)
println('')
println('Currently using the Cecret workflow for use with amplicon Illumina library prep on MiSeq with a corresponding reference genome.\n')
println('Author: Erin Young')
println('email: eriny@utah.gov')
println('')

println('Cecret is named after a real lake!')
println('Visit https://www.alltrails.com/trail/us/utah/cecret-lake-trail to learn more.')
println('Not everyone can visit in person, so here is some ASCII art of nucleotides in lake forming a consensus sequence.')
println('                    _________ ______')
println('               _ /      G    A   T   \\_____')
println('          __/    C      C A    G      T  C \\')
println('        /    G     A   T   T  A   G  G    T  \\_')
println('        | G       G  C   A            G   T     \\')  
println('        \\      A     C     G   A   T    A  G  T  \\__')
println('         \\_           C       G    ____ _____ __ C  \\________')
println('            \\__T______ ___________/                \\ C T G A G G T C G A T A') 
println('')
println('')

//# copying the confit template and ending the workflow
params.config_file                          = false
if (params.config_file) {
  def src = new File("${workflow.projectDir}/configs/cecret_config_template.config")
  def dst = new File("${workflow.launchDir}/edit_me.config")
  dst << src.text
  println('A config file can be found at ' + workflow.launchDir + '/edit_me.config')
  exit 0
}

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

//# Starting the workflow --------------------------------------------------------------

nextflow.enable.dsl = 2

//# params and their default values

//# params for inputs including fasta and fastq files
params.sample_sheet                         = ''
params.reads                                = workflow.launchDir + '/reads'
params.single_reads                         = workflow.launchDir + '/single_reads'
params.fastas                               = workflow.launchDir + '/fastas'
params.multifastas                          = workflow.launchDir + '/multifastas'
params.nanopore                             = workflow.launchDir + '/nanopore'
params.sra_accessions                       = []
params.primer_bed                           = ''
params.amplicon_bed                         = ''
params.primer_set                           = 'ncov_V4'
params.reference_genome                     = ''
params.gff                                  = ''

//# input checks ---------------------------------------------------------------------

//# Ensuring that reads and single_reads are not set to the same directory
if ( params.reads == params.single_reads ) {
  println('\'params.reads\' and \'params.single_reads\' cannot point to the same directory!')
  println('\'params.reads\' is set to ' + params.reads)
  println('\'params.single_reads\' is set to ' + params.single_reads)
  exit 1
}

if ( params.fastas == params.multifastas ) {
  println('\'params.fastas\' and \'params.multifastas\' cannot point to the same directory!')
  println('\'params.fastas\' is set to ' + params.fastas)
  println('\'params.multifastas\' is set to ' + params.multifastas)
  exit 1
}

//# run params ---------------------------------------------------------------------

//# outdir params
params.outdir                               = workflow.launchDir + '/cecret'
params.species                              = 'sarscov2'
println('The files and directory for results is ' + params.outdir)
println('The species used to determine default variables and subworkflows is ' + params.species)

//# roughly grouping cpu usage
params.maxcpus                              = 8
params.medcpus                              = 4
println('The maximum number of CPUS used in this workflow is ' + params.maxcpus)

//# specifying the core workflow
params.trimmer                              = 'ivar'
params.cleaner                              = 'seqyclean'
params.aligner                              = 'bwa'
params.msa                                  = 'mafft'

//# to toggle off processes
params.aci                                  = true
params.bcftools_variants                    = true
params.fastqc                               = true
params.ivar_variants                        = true
params.samtools_stats                       = true
params.samtools_coverage                    = true
params.samtools_depth                       = true
params.samtools_flagstat                    = true
params.samtools_ampliconstats               = true
params.samtools_plot_ampliconstats          = true
params.markdup                              = false
params.filter                               = false
params.multiqc                              = true

//# for optional route of tree generation and counting snps between samples
params.relatedness                          = false
params.snpdists                             = true
params.iqtree2                              = true

//# parameters for processes with their default values
params.artic_options                        = '--normalise 200 --skip-nanopolish --medaka --medaka-model r941_min_high_g360'
params.artic_read_filtering_options         = '--min-length 400 --max-length 700'
params.bcftools_variants_options            = ''
params.fastp_options                        = ''
params.fastqc_options                       = ''
params.filter_options                       = ''
params.iqtree2_options                      = '-ninit 2 -n 2 -me 0.05 -m GTR'
params.ivar_consensus_options               = '-q 20 -t 0.6 -n N'
params.ivar_trim_options                    = ''
params.ivar_variants_options                = '-q 20 -t 0.6'
params.minimap2_options                     = '-K 20M'
params.minimum_depth                        = 100
params.mpileup_depth                        = 8000
params.multiqc_options                      = ''
params.mafft_options                        = '--maxambiguous 0.5'
params.samtools_ampliconclip_options        = ''
params.samtools_coverage_options            = ''
params.samtools_flagstat_options            = ''
params.samtools_depth_options               = ''
params.samtools_stats_options               = ''
params.samtools_ampliconstats_options       = '--max-amplicon-length 3000 --max-amplicons 3000'
params.samtools_plot_ampliconstats_options  = '-size 1200,900 -size2 1200,900 -size3 1200,900'
params.samtools_markdup_options             = ''
params.samtools_fixmate_options             = ''
params.seqyclean_contaminant_file           = '/Adapters_plus_PhiX_174.fasta'
params.seqyclean_options                    = '-minlen 25 -qual'
params.snpdists_options                     = '-c'

//# for optional contamination determination
params.kraken2                              = false
params.kraken2_db                           = false
params.kraken2_options                      = ''

//# for using an included version of nextclade dataset
params.download_nextclade_dataset           = true

//# organism specific
params.freyja                               = true
params.freyja_aggregate                     = true
params.nextclade                            = true
params.pangolin                             = true
params.vadr                                 = true

params.pangolin_options                     = ''
params.vadr_mdir                            = '/opt/vadr/vadr-models'
params.nextclade_options                    = ''
params.nextalign_options                    = '--include-reference'
params.freyja_variants_options              = ''
params.freyja_demix_options                 = ''
params.freyja_boot_options                  = '--nb 1000'
params.freyja_aggregate_options             = ''
params.freyja_plot_options                  = ''
params.freyja_plot_filetype                 = 'png'

//# Specifying some species-specific params
if ( params.species == 'sarscov2' ) {
  params.nextclade_dataset                  = 'sars-cov-2'
  params.vadr_options                       = '--split --glsearch -s -r --nomisc --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn'
  params.vadr_reference                     = 'sarscov2'
  params.vadr_trim_options                  = '--minlen 50 --maxlen 30000'
  params.iqtree2_outgroup                   = 'MN908947'
} else if ( params.species == 'mpx' ) {
  params.nextclade_dataset                  = 'hMPXV'
  params.vadr_options                       = '--split --glsearch -s -r --nomisc --r_lowsimok --r_lowsimxd 100 --r_lowsimxl 2000 --alt_pass discontn,dupregin'
  params.vadr_reference                     = 'mpxv'
  params.vadr_trim_options                  = '--minlen 50 --maxlen 210000'
  params.iqtree2_outgroup                   = 'NC_063383'
} else {
  params.nextclade_dataset                  = ''
  params.vadr_options                       = ''
  params.vadr_reference                     = ''
  params.vadr_trim_options                  = ''
  params.iqtree2_outgroup                   = ''
}

//# Adding in subworkflows
include { fasta_prep ; summary } from './modules/cecret.nf'      addParams(params)
include { cecret }               from './subworkflows/cecret.nf' addParams(params)
include { qc }                   from './subworkflows/qc'        addParams(params)
include { msa }                  from './subworkflows/msa'       addParams(params)
include { multiqc_combine }      from './modules/multiqc'        addParams(params)
include { mpx }                  from './subworkflows/mpx'       addParams(params)                                 
include { mpx as other }         from './subworkflows/mpx'       addParams(params)
include { sarscov2 }             from './subworkflows/sarscov2'  addParams(params)
include { test }                 from './subworkflows/test'      addParams(params) 

//# Now that everything is defined (phew!), the workflow can begin ---------------------------------------------------
//# getting input files
if ( params.sample_sheet ) { 
  Channel
    .fromPath(params.sample_sheet, type: 'file')
    .view { "Sample sheet found : ${it}" }
    .splitCsv( header: true, sep: ',' )
    .map { row -> tuple( "${row.sample}", [ file("${row.fastq_1}"), file("${row.fastq_2}") ]) }
    .branch {
      single     : it[1] =~ /single/
      multifasta : it[1] =~ /multifasta/
      fasta      : it[1] =~ /fasta/
      ont        : it[1] =~ /ont/       
      paired     : true 
    }
    .set { inputs }
  
  ch_paired_reads = inputs.paired.map{     it -> tuple(it[0], it[1],    'paired')}
  ch_single_reads = inputs.single.map{     it -> tuple(it[0], it[1][0], 'single')}
  ch_fastas       = inputs.fasta.map{      it -> tuple(it[0], it[1])}
  ch_multifastas  = inputs.multifasta.map{ it -> tuple(it[1])}
  ch_nanopore     = inputs.ont.map{        it -> tuple(it[0], it[1])}
} else {
  Channel
    .fromFilePairs(["${params.reads}/*_R{1,2}*.{fastq,fastq.gz,fq,fq.gz}",
                    "${params.reads}/*{1,2}*.{fastq,fastq.gz,fq,fq.gz}"], size: 2 )
    .unique()
    .map { reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1], 'paired' ) }
    .set { ch_paired_reads }

  Channel
    .fromPath("${params.single_reads}/*.{fastq,fastq.gz,fq,fq.gz}")
    .map { reads -> tuple(reads.simpleName, reads, 'single' ) }
    .set { ch_single_reads }

    Channel
    .fromPath("${params.nanopore}/*.{fastq,fastq.gz,fq,fq.gz}")
    .map { reads -> tuple(reads.simpleName, reads ) }
    .set { ch_nanopore }

  Channel
    .fromPath("${params.fastas}/*{.fa,.fasta,.fna}", type:'file')
    .map { fasta -> tuple(fasta.baseName, fasta ) }
    .set { ch_fastas }

  Channel
    .fromPath("${params.multifastas}/*{.fa,.fasta,.fna}", type:'file')
    .set { ch_multifastas }
}

ch_sra_accessions = Channel.from( params.sra_accessions )

//# Checking for input files and giving an explanatory message if none are found
ch_paired_reads
  .mix(ch_single_reads)
  .mix(ch_fastas)
  .mix(ch_multifastas)
  .mix(ch_sra_accessions)
  .mix(ch_nanopore)
  .ifEmpty{
    println('FATAL : No input files were found!')
    println("No paired-end fastq files were found at ${params.reads}. Set 'params.reads' to directory with paired-end reads")
    println("No single-end fastq files were found at ${params.single_reads}. Set 'params.single_reads' to directory with single-end reads")
    println("No fasta files were found at ${params.fastas}. Set 'params.fastas' to directory with fastas.")
    println("No multifasta files were found at ${params.multifastas}. Set 'params.multifastas' to directory with multifastas.")
    println("No sample sheet was fount at ${params.sample_sheet}. Set 'params.sample_sheet' to sample sheet file.")
    exit 1
}

//# getting a reference genome file
if (params.reference_genome){
  Channel
    .fromPath(params.reference_genome, type:'file')
    .ifEmpty{
      println("No reference genome was selected. Set with 'params.reference_genome'")
      exit 1
    }
    .set { ch_reference_genome }
} else {
  if ( params.species == 'sarscov2' ) {
    ch_reference_genome = Channel.fromPath(workflow.projectDir + '/genomes/MN908947.3.fasta', type: 'file')
  } else if ( params.species == 'mpx') {
    ch_reference_genome = Channel.fromPath(workflow.projectDir + '/genomes/NC_063383.1.fasta', type: 'file')
  } else {
    println("No reference genome was selected. Set with 'params.reference_genome'")
    println("Or set species to one with an included genome ('sarscov2' or 'mpx')")
    exit 1
    ch_reference_genome = Channel.empty()
  } 
}
ch_reference_genome.view { "Reference Genome : $it"}

//# getting the gff file for ivar variants
if ( params.ivar_variants ) {
  if (params.gff) {
    Channel
      .fromPath(params.gff, type:'file')
      .ifEmpty{
        println("No gff file was selected. Set with 'params.reference_genome'")
        exit 1
      }
      .set { ch_gff_file }
  } else {
    if ( params.species == 'sarscov2' ) {
      ch_gff_file = Channel.fromPath(workflow.projectDir + '/genomes/MN908947.3.gff', type: 'file')
    } else if ( params.species == 'mpx') {
      ch_gff_file = Channel.fromPath(workflow.projectDir + '/genomes/NC_063383.1.gff3', type: 'file')
    } else {
      println("No gff file was selected. Set with 'params.gff'")
      println("Or set 'params.species' to one with an included genome ('sarscov2' or 'mpx')")
      println("Or bypass this message completely by setting 'params.ivar_variants = False'")
      exit 1
      ch_gff_file = Channel.empty()
    } 
  }
}
ch_gff_file.view { "GFF file : $it"}

//# channels of included files
included_primers     = [
  workflow.projectDir + '/schema/midnight_idt_V1_SARS-CoV-2.primer.bed',
  workflow.projectDir + '/schema/midnight_ont_V1_SARS-CoV-2.primer.bed',
  workflow.projectDir + '/schema/midnight_ont_V2_SARS-CoV-2.primer.bed',
  workflow.projectDir + '/schema/midnight_ont_V3_SARS-CoV-2.primer.bed',
  workflow.projectDir + '/schema/ncov_V3_nCoV-2019.primer.bed',
  workflow.projectDir + '/schema/ncov_V4_SARS-CoV-2.primer.bed',
  workflow.projectDir + '/schema/ncov_V4.1_SARS-CoV-2.primer.bed',
  workflow.projectDir + '/schema/ncov_V5.3.2_SARS-CoV-2.primer.bed',
  workflow.projectDir + '/schema/mpx_idt_primer.bed',
  workflow.projectDir + '/schema/mpx_primalseq_primer.bed'
  ]
included_amplicons = [
  workflow.projectDir + '/schema/midnight_idt_V1_SARS-CoV-2.insert.bed',
  workflow.projectDir + '/schema/midnight_ont_V1_SARS-CoV-2.insert.bed',
  workflow.projectDir + '/schema/midnight_ont_V2_SARS-CoV-2.insert.bed',
  workflow.projectDir + '/schema/midnight_ont_V3_SARS-CoV-2.insert.bed',
  workflow.projectDir + '/schema/ncov_V3_nCoV-2019.insert.bed',
  workflow.projectDir + '/schema/ncov_V4_SARS-CoV-2.insert.bed',
  workflow.projectDir + '/schema/ncov_V4.1_SARS-CoV-2.insert.bed',
  workflow.projectDir + '/schema/ncov_V5.3.2_SARS-CoV-2.insert.bed',
  workflow.projectDir + '/schema/mpx_idt_insert.bed',
  workflow.projectDir + '/schema/mpx_primalseq_insert.bed'
]

ch_primers                = Channel.fromPath(included_primers,   type: 'file')
ch_amplicons              = Channel.fromPath(included_amplicons, type: 'file')

available_primer_sets = [
  'midnight_idt_V1', 
  'midnight_ont_V1', 
  'midnight_ont_V2', 
  'midnight_ont_V3', 
  'ncov_V3', 
  'ncov_V4', 
  'ncov_V4.1', 
  'ncov_V5.3.2', 
  'mpx_primalseq', 
  'mpx_idt'
  ]

if ( params.trimmer != 'none' ) {
  //# Getting the primer file
  if (params.primer_bed) {
    Channel
      .fromPath(params.primer_bed, type:'file')
      .ifEmpty{
        println("A bedfile for primers is required. Set with 'params.primer_bed'.")
        exit 1
      }
      .set { ch_primer_bed } 
  } else if ( params.primer_set in available_primer_sets ) {
    Channel
      .fromPath( included_primers )
      .branch{ 
        match : it =~ /${params.primer_set}_*/
        }
      .first()
      .set { ch_primer_bed } 
  } else {
    println("No primers were found!")
    println("Set primer schema with 'params.primer_bed' or specify to 'none' if primers were not used")
    println("Or use included primer set by setting 'params.primer_set' to one of $available_primer_sets")
    exit 1
    ch_primer_bed = Channel.empty()
  }
  ch_primer_bed.view { "Primer BedFile : $it"}

  //# Getting the amplicon bedfile
  if ( params.aci ) {
    if (params.amplicon_bed) {
      Channel
        .fromPath(params.amplicon_bed, type:'file')
        .ifEmpty{
          println("A bedfile for amplicons is required. Set with 'params.amplicon_bed'.")
          println("Or set params.aci = False to skip this.")
          exit 1
        }
        .set { ch_amplicon_bed } 
    } else if ( params.primer_set in available_primer_sets ) {
      Channel
        .fromPath( included_amplicons )
        .branch{ 
          match : it =~ /${params.primer_set}_*/
          }
        .first()
        .set { ch_amplicon_bed } 
    } else {
      println("An amplicon bedfile wasn't found!")
      println("Set amplicon schema with 'params.amplicon_bed'")
      println("Or use included primer set by setting 'params.primer_set' to one of $available_primer_sets")
      println("Or set params.aci = False to skip this.")
      exit 1
      ch_amplicon_bed = Channel.empty()
    }
    ch_amplicon_bed.view { "Amplicon BedFile : $it"}
  } else {
    ch_amplicon_bed = Channel.empty()
  }
} else {
  ch_primer_bed   = Channel.empty()
  ch_amplicon_bed = Channel.empty()
}

//# scripts for legacy reasons
ch_combine_results_script = Channel.fromPath("${workflow.projectDir}/bin/combine_results.py",  type:'file')
ch_freyja_script          = Channel.fromPath("${workflow.projectDir}/bin/freyja_graphs.py",    type:'file')

if ( params.kraken2_db ) {
  Channel
    .fromPath(params.kraken2_db, type:'dir')
    .view { "Kraken2 database : $it" }
    .set{ ch_kraken2_db }
} else {
  ch_kraken2_db = Channel.empty()
}

if ( ! params.download_nextclade_dataset ) {
  Channel
    .fromPath(params.predownloaded_nextclade_dataset)
    .ifEmpty{
      println("Dataset file could not be found at ${params.predownloaded_nextclade_dataset}.")
      println("Please set nextclade dataset file with 'params.predownloaded_nextclade_dataset'")
      exit 1
    }
    .set { ch_nextclade_dataset }
} else {
  ch_nextclade_dataset = Channel.empty()
}

ch_paired_reads
  .mix(ch_single_reads)
  .unique()
  .set { ch_reads }

ch_paired_reads.view { "Paired-end Fastq files found : ${it[0]}" }
ch_single_reads.view { "Fastq files found : ${it[0]}" }
ch_fastas.view       { "Fasta file found : ${it[0]}" }
ch_multifastas.view  { "MultiFasta file found : ${it}" }
ch_reads.ifEmpty     { println("No fastq or fastq.gz files were found at ${params.reads} or ${params.single_reads}") }

workflow CECRET {
    ch_for_dataset = Channel.empty()
    ch_for_version = Channel.from("Cecret version", workflow.manifest.version).collect()

    if ( ! params.sra_accessions.isEmpty() ) { 
      test(ch_sra_accessions)
      ch_reads = ch_reads.mix(test.out.reads)
    } 

    fasta_prep(ch_fastas)

    cecret(ch_reads, ch_nanopore, ch_reference_genome, ch_primer_bed)

    qc(ch_reads,
      cecret.out.clean_reads,
      ch_kraken2_db,
      cecret.out.sam,
      cecret.out.trim_bam,
      ch_reference_genome,
      ch_gff_file,
      ch_amplicon_bed,
      ch_primer_bed)

    ch_for_multiqc = cecret.out.for_multiqc.mix(qc.out.for_multiqc)
    ch_for_summary = qc.out.for_summary

    if ( params.species == 'sarscov2' ) {
      sarscov2(fasta_prep.out.fastas.mix(ch_multifastas).mix(cecret.out.consensus), cecret.out.trim_bam, ch_reference_genome, ch_nextclade_dataset, ch_freyja_script)
      
      ch_for_multiqc = ch_for_multiqc.mix(sarscov2.out.for_multiqc)
      ch_for_dataset = sarscov2.out.dataset
      ch_for_summary = ch_for_summary.mix(sarscov2.out.for_summary)
    
    } else if ( params.species == 'mpx') {
      mpx(fasta_prep.out.fastas.mix(ch_multifastas).mix(cecret.out.consensus), ch_nextclade_dataset)
      
      ch_for_multiqc = ch_for_multiqc.mix(mpx.out.for_multiqc)
      ch_for_dataset = mpx.out.dataset
      ch_for_summary = ch_for_summary.mix(mpx.out.for_summary)

    } else if ( params.species == 'other') {
      other(fasta_prep.out.fastas.concat(ch_multifastas).mix(cecret.out.consensus), ch_nextclade_dataset)
      
      ch_for_multiqc = ch_for_multiqc.mix(other.out.for_multiqc)
      ch_for_dataset = other.out.dataset
      ch_for_summary = ch_for_summary.mix(other.out.for_summary)

    } 

    if ( params.relatedness ) { 
      msa(fasta_prep.out.fastas.concat(ch_multifastas).concat(cecret.out.consensus), ch_reference_genome, ch_for_dataset) 

      tree      = msa.out.tree
      alignment = msa.out.msa
      matrix    = msa.out.matrix

    } else {
      tree      = Channel.empty()
      alignment = Channel.empty()
      matrix    = Channel.empty()
    }

    multiqc_combine(ch_for_multiqc.collect())

    summary(
      ch_for_summary.mix(fasta_prep.out.fastas).mix(cecret.out.consensus).collect().map{it -> tuple([it])}
        .combine(ch_combine_results_script)
        .combine(ch_for_version.mix(cecret.out.for_version).collect().map{it -> tuple([it])})
        .combine(multiqc_combine.out.files.ifEmpty([]).map{it -> tuple([it])}))

  emit:
    bam       = cecret.out.trim_bam
    consensus = fasta_prep.out.fastas.mix(ch_multifastas).mix(cecret.out.consensus).collect()
    tree      = tree
    alignment = alignment
    matrix    = matrix
}

workflow {
  CECRET ()
}

workflow.onComplete {
  println("Pipeline completed at: $workflow.complete")
  println("A summary of results can be found in a comma-delimited file: ${params.outdir}/cecret_results.csv")
  println("A summary of results can be found in a tab-delimited file: ${params.outdir}/cecret_results.txt")
  println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
}
