#!/usr/bin/env nextflow

println('Currently using the Cecret workflow for use with amplicon Illumina library prep on MiSeq with a corresponding reference genome.\n')
println('Author: Erin Young')
println('email: eriny@utah.gov')
println("Version: ${workflow.manifest.version}")
println('')

params.config_file                          = false
if (params.config_file) {
  def src = new File("${workflow.projectDir}/configs/cecret_config_template.config")
  def dst = new File("${workflow.launchDir}/edit_me.config")
  dst << src.text
  println("A config file can be found at ${workflow.launchDir}/edit_me.config")
  exit 0
}

nextflow.enable.dsl = 2

//# params

//# input params
params.reads                                = workflow.launchDir + '/reads'
params.single_reads                         = workflow.launchDir + '/single_reads'
params.fastas                               = workflow.launchDir + '/fastas'
params.multifastas                          = workflow.launchDir + '/multifastas'

if ( params.reads == params.single_reads ) {
  println("'params.reads' and 'params.single_reads' cannot point to the same directory!")
  println("'params.reads' is set to " + params.reads)
  println("'params.single_reads' is set to " + params.single_reads)
  exit 1
}

if ( params.fastas == params.multifastas ) {
  println("'params.fastas' and 'params.multifastas' cannot point to the same directory!")
  println("'params.fastas' is set to " + params.fastas)
  println("'params.multifastas' is set to " + params.multifastas)
  exit 1
}

//# outdir params
params.outdir                               = workflow.launchDir + '/cecret'

params.maxcpus                              = 8
params.medcpus                              = 4
println("The maximum number of CPUS used in this workflow is ${params.maxcpus}")

//# default reference files for SARS-CoV-2 or MPX (part of the github repository)
params.species                              = 'sarscov2'
if (params.species        == 'sarscov2' ) {
  params.reference_genome                   = workflow.projectDir + '/configs/MN908947.3.fasta'
  params.gff                                = workflow.projectDir + '/configs/MN908947.3.gff'
  println("Using the subworkflow for SARS-CoV-2")
} else if (params.species == 'mpx') {
  params.reference_genome                   = workflow.projectDir + '/configs/NC_063383.1.fasta'
  params.gff                                = workflow.projectDir + '/configs/NC_063383.1.gff3'
  println("Using the subworkflow for Monkeypox Virus")
} else {
  params.reference_genome                   = ''
  params.gff                                = ''
}

params.primer_set                           = 'ncov_V4'
if ( params.primer_set        == 'ncov_V3' ) {
  params.primer_bed                         = workflow.projectDir + '/configs/artic_V3_nCoV-2019.primer.bed'
  params.amplicon_bed                       = workflow.projectDir + '/configs/artic_V3_nCoV-2019.insert.bed'
} else if ( params.primer_set == 'ncov_V4' ) {
  params.primer_bed                         = workflow.projectDir + '/configs/artic_V4_SARS-CoV-2.primer.bed'
  params.amplicon_bed                       = workflow.projectDir + '/configs/artic_V4_SARS-CoV-2.insert.bed'
} else if ( params.primer_set == 'ncov_V4.1' ) {
  params.primer_bed                         = workflow.projectDir + '/configs/artic_V4.1_SARS-CoV-2.primer.bed'
  params.amplicon_bed                       = workflow.projectDir + '/configs/artic_V4.1_SARS-CoV-2.insert.bed'
} else if ( params.primer_set == 'mpx_idt' ) {
  params.primer_bed                         = workflow.projectDir + '/configs/mpx_idt_primer.bed'
  params.amplicon_bed                       = workflow.projectDir + '/configs/mpx_idt_insert.bed'
} else {
  println("!{params.primer_set} has not been defined as an acceptable value for 'params.primer_set'.")
  println('Current acceptable values are' )
  println("SARS-CoV-2 artic primer V3 : 'params.primer_set' = 'ncov_V3'" )
  println("SARS-CoV-2 artic primer V4 : 'params.primer_set' = 'ncov_V4'" )
  println("SARS-CoV-2 artic primer V4.1 (Version 4 with spike in) : 'params.primer_set' = 'ncov_V4.1'" )
  println("Monkeypox IDT primer : 'params.primer_set' = 'mpx_idt'" )
  exit 1
}

//# specifying the core workflow
params.trimmer                              = 'ivar'
params.cleaner                              = 'seqyclean'
params.aligner                              = 'bwa'
params.msa                                  = 'mafft'

//# to toggle off processes
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
params.bedtools_multicov                    = true
params.kraken2                              = false
params.filter                               = false
params.multiqc                              = true

//# for optional route of tree generation and counting snps between samples
params.relatedness                          = false
params.snpdists                             = true
params.iqtree2                              = true

//# parameters for processes
params.fastqc_options                       = ''
params.seqyclean_contaminant_file           = '/Adapters_plus_PhiX_174.fasta'
params.seqyclean_options                    = '-minlen 25 -qual'
params.fastp_options                        = ''
params.minimap2_options                     = '-K 20M'
params.filter_options                       = ''
params.ivar_trim_options                    = ''
params.samtools_ampliconclip_options        = ''
params.minimum_depth                        = 100
params.mpileup_depth                        = 8000
params.ivar_variants_options                = '-q 20 -t 0.6'
params.ivar_consensus_options               = '-q 20 -t 0.6 -n N'
params.kraken2_options                      = ''
params.bedtools_multicov_options            = '-f .1'
params.bcftools_variants_options            = ''
params.samtools_coverage_options            = ''
params.samtools_flagstat_options            = ''
params.samtools_depth_options               = ''
params.samtools_stats_options               = ''
params.samtools_ampliconstats_options       = ''
params.samtools_plot_ampliconstats_options  = '-size 1200,900 -size2 1200,900 -size3 1200,900'
params.samtools_markdup_options             = ''
params.samtools_fixmate_options             = ''
params.mafft_options                        = '--maxambiguous 0.5'
params.snpdists_options                     = '-c'
params.iqtree2_options                      = '-ninit 2 -n 2 -me 0.05 -m GTR'
params.multiqc_options                      = ''

//# for optional contamination determination
params.kraken2_db                           = false

//# organism specific
params.nextclade                            = true
params.pangolin                             = true
params.vadr                                 = true
params.freyja                               = true
params.freyja_aggregate                     = true

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

if ( params.species == 'sarscov2' ) {
  params.nextclade_dataset                  = 'sars-cov-2'
  params.vadr_options                       = '--split --glsearch -s -r --nomisc --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn'
  params.vadr_reference                     = 'sarscov2'
  params.vadr_trim_options                  = '--minlen 50 --maxlen 30000'
  params.kraken2_organism                   = 'Severe acute respiratory syndrome-related coronavirus'
  params.iqtree2_outgroup                   = 'MN908947'
} else if ( params.species == 'mpx' ) {
  params.nextclade_dataset                  = 'hMPXV'
  params.vadr_options                       = '--split --glsearch -s -r --nomisc --r_lowsimok --r_lowsimxd 100 --r_lowsimxl 2000 --alt_pass discontn,dupregin'
  params.vadr_reference                     = 'mpxv'
  params.vadr_trim_options                  = '--minlen 50 --maxlen 210000'
  params.kraken2_organism                   = 'Monkeypox virus'
  params.iqtree2_outgroup                   = 'NC_063383'
} else {
  params.nextclade_dataset                  = ''
  params.vadr_options                       = ''
  params.vadr_reference                     = ''
  params.vadr_trim_options                  = ''
  params.kraken2_organis                    = '.'
}

include { fasta_prep ; summary; combine_results } from './modules/cecret.nf'      addParams(fastqc: params.fastqc,
                                                                                            trimmer: params.trimmer,
                                                                                            cleaner: params.cleaner,
                                                                                            samtools_coverage: params.samtools_coverage,
                                                                                            samtools_depth: params.samtools_depth,
                                                                                            minimum_depth: params.minimum_depth,
                                                                                            samtools_stats: params.samtools_stats,
                                                                                            samtools_ampliconstats: params.samtools_ampliconstats,
                                                                                            kraken2: params.kraken2,
                                                                                            ivar_variants: params.ivar_variants,
                                                                                            bcftools_variants: params.bcftools_variants,
                                                                                            bedtools_multicov: params.bedtools_multicov)
include { cecret }                                from './subworkflows/cecret.nf' addParams(cleaner: params.cleaner,
                                                                                            aligner: params.aligner,
                                                                                            trimmer: params.trimmer,
                                                                                            seqyclean_options: params.seqyclean_options,
                                                                                            seqyclean_contaminant_file: params.seqyclean_contaminant_file,
                                                                                            fastp_options: params.fastp_options,
                                                                                            minimap2_options: params.minimap2_options,
                                                                                            ivar_trim_options: params.ivar_trim_options,
                                                                                            samtools_ampliconclip_options: params.samtools_ampliconclip_options,
                                                                                            minimum_depth: params.minimum_depth,
                                                                                            mpileup_depth: params.mpileup_depth,
                                                                                            samtools_fixmate_options: params.samtools_fixmate_options,
                                                                                            samtools_markdup_options: params.samtools_markdup_options,
                                                                                            ivar_consensus_options: params.ivar_consensus_options)
include { qc }                                    from './subworkflows/qc'        addParams(trimmer: params.trimmer,
                                                                                            fastqc: params.fastqc,
                                                                                            fastqc_options: params.fastqc_options,
                                                                                            kraken2: params.kraken2,
                                                                                            kraken2_options: params.kraken2_options,
                                                                                            kraken2_organism: params.kraken2_organism,
                                                                                            bcftools_variants: params.bcftools_variants,
                                                                                            bcftools_variants_options: params.bcftools_variants_options,
                                                                                            ivar_variants: params.ivar_variants,
                                                                                            ivar_variants_options: params.ivar_variants_options,
                                                                                            bedtools_multicov: params.bedtools_multicov,
                                                                                            bedtools_multicov_options: params.bedtools_multicov_options,
                                                                                            samtools_stats: params.samtools_stats,
                                                                                            samtools_stats_options: params.samtools_stats_options,
                                                                                            samtools_coverage: params.samtools_coverage,
                                                                                            samtools_coverage_options: params.samtools_coverage_options,
                                                                                            samtools_flagstat: params.samtools_flagstat,
                                                                                            samtools_flagstat_options: params.samtools_flagstat_options,
                                                                                            samtools_ampliconstats: params.samtools_ampliconstats,
                                                                                            samtools_ampliconstats_options: params.samtools_ampliconstats_options,
                                                                                            samtools_plot_ampliconstats: params.samtools_plot_ampliconstats,
                                                                                            samtools_plot_ampliconstats_options: params.samtools_plot_ampliconstats_options)
include { msa }                                   from './subworkflows/msa'       addParams(msa: params.msa,
                                                                                            nextalign_options: params.nextalign_options,
                                                                                            mafft_options: params.mafft_options,
                                                                                            iqtree2: params.iqtree2,
                                                                                            iqtree2_options: params.iqtree2_options,
                                                                                            iqtree2_outgroup: params.iqtree2_outgroup,
                                                                                            snpdists: params.snpdists,
                                                                                            snpdists_options: params.snpdists_options)
include { multiqc_combine }                       from './modules/multiqc'        addParams(multiqc: params.multiqc,
                                                                                            multiqc_options: params.multiqc_options)
include { mpx }                                   from './subworkflows/mpx'       addParams(vadr: params.vadr,
                                                                                            vadr_options: params.vadr_options,
                                                                                            vadr_reference: params.vadr_reference,
                                                                                            vadr_mdir: params.vadr_mdir,
                                                                                            nextclade: params.nextclade,
                                                                                            nextclade_options: params.nextclade_options,
                                                                                            nextclade_dataset: params.nextclade_dataset)                                 
include { mpx as other }                          from './subworkflows/mpx'       addParams(vadr: params.vadr,
                                                                                            vadr_options: params.vadr_options,
                                                                                            vadr_reference: params.vadr_reference,
                                                                                            vadr_mdir: params.vadr_mdir,
                                                                                            nextclade: params.nextclade,
                                                                                            nextclade_options: params.nextclade_options,
                                                                                            nextclade_dataset: params.nextclade_dataset)
include { sarscov2 }                              from './subworkflows/sarscov2'  addParams(vadr: params.vadr,
                                                                                            vadr_options: params.vadr_options,
                                                                                            vadr_reference: params.vadr_reference,
                                                                                            vadr_mdir: params.vadr_mdir,
                                                                                            pangolin: params.pangolin,
                                                                                            pangolin_options: params.pangolin_options,
                                                                                            nextclade: params.nextclade,
                                                                                            nextclade_options: params.nextclade_options,
                                                                                            nextclade_dataset: params.nextclade_dataset,
                                                                                            freyja: params.freyja,
                                                                                            freyja_variants_options: params.freyja_variants_options,
                                                                                            freyja_demix_options: params.freyja_demix_options,
                                                                                            freyja_aggregate: params.freyja_aggregate,
                                                                                            freyja_aggregate_options: params.freyja_aggregate_options,
                                                                                            freyja_plot_options: params.freyja_plot_options)

//# getting input files
Channel
  .fromFilePairs(["${params.reads}/*_R{1,2}*.{fastq,fastq.gz,fq,fq.gz}",
                  "${params.reads}/*{1,2}*.{fastq,fastq.gz,fq,fq.gz}"], size: 2 )
  .unique()
  .map { reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1], "paired" ) }
  .set { paired_reads }

Channel
  .fromPath("${params.single_reads}/*.{fastq,fastq.gz,fq,fq.gz}")
  .map { reads -> tuple(reads.simpleName, reads, "single" ) }
  .set { single_reads }

Channel
  .fromPath("${params.fastas}/*{.fa,.fasta,.fna}", type:'file')
  .map { fasta -> tuple(fasta.baseName, fasta ) }
  .set { fastas }

multifastas = Channel.fromPath("${params.multifastas}/*{.fa,.fasta,.fna}", type:'file')

//# Checking for input files and giving an explanatory message if none are found
paired_reads
  .mix(single_reads)
  .mix(fastas)
  .mix(multifastas)
  .ifEmpty{
    println('FATAL : No input files were found!')
    println("No paired-end fastq files were found at ${params.reads}. Set 'params.reads' to directory with paired-end reads")
    println("No single-end fastq files were found at ${params.single_reads}. Set 'params.single_reads' to directory with single-end reads")
    println("No fasta files were found at ${params.fastas}. Set 'params.fastas' to directory with fastas.")
    println("No multifasta files were found at ${params.multifastas}. Set 'params.multifastas' to directory with multifastas.")
    exit 1
}

//# getting reference files
Channel
  .fromPath(params.reference_genome, type:'file')
  .ifEmpty{
    println("No reference genome was selected. Set with 'params.reference_genome'")
    exit 1
  }
  .view { "Reference Genome : $it"}
  .set { reference_genome }

gff_file = params.ivar_variants
  ? Channel.fromPath(params.gff, type:'file').view { "GFF file for Reference Genome : $it"}
  : Channel.empty()

if ( params.trimmer != 'none' ) {
  Channel
    .fromPath(params.primer_bed, type:'file')
    .ifEmpty{
      println("A bedfile for primers is required. Set with 'params.primer_bed'.")
      exit 1
    }
    .view { "Primer BedFile : $it"}
    .set { primer_bed }

  amplicon_bed = params.bedtools_multicov
    ? Channel.fromPath(params.amplicon_bed, type:'file').view { "Amplicon BedFile : $it"}
    : Channel.empty()
} else {
  primer_bed = Channel.empty()
  amplicon_bed = Channel.empty()
}

kraken2_db = params.kraken2_db
  ? Channel.fromPath(params.kraken2_db, type:'dir').view { "Kraken2 database : $it" }
  : Channel.empty()

//# getting scripts
combine_results_script = Channel.fromPath("${workflow.projectDir}/bin/combine_results.py", type:'file')

// This is where the results will be
println('The files and directory for results is ' + params.outdir)
println("A table summarizing results will be created: ${params.outdir}/cecret_results.csv\n")

paired_reads
  .mix(single_reads)
  .unique()
  .set { reads }

workflow {
  paired_reads.view { "Paired-end Fastq files found : ${it[0]}" }
  single_reads.view { "Fastq files found : ${it[0]}" }
  fastas.view { "Fasta file found : ${it[0]}" }
  multifastas.view { "MultiFasta file found : ${it}" }
  reads.ifEmpty{ println("No fastq or fastq.gz files were found at ${params.reads} or ${params.single_reads}") }

  //  combine_results_script
  fasta_prep(fastas)

  cecret(reads,reference_genome,primer_bed)
  qc(reads,
    cecret.out.clean_type,
    kraken2_db,
    cecret.out.sam,
    cecret.out.bam,
    cecret.out.bam_bai,
    reference_genome,
    gff_file,
    amplicon_bed,
    primer_bed)

  if ( params.species == 'sarscov2' ) {
    sarscov2(fasta_prep.out.fastas.mix(multifastas).mix(cecret.out.consensus), cecret.out.bam, reference_genome)
    pangolin_file   = sarscov2.out.pangolin_file
    nextclade_file  = sarscov2.out.nextclade_file
    vadr_file       = sarscov2.out.vadr_file
    freyja_file     = sarscov2.out.freyja_file
    dataset         = sarscov2.out.dataset 
  } else if ( params.species == 'mpx') {
    mpx(fasta_prep.out.fastas.mix(multifastas).mix(cecret.out.consensus))
    pangolin_file   = Channel.empty()
    freyja_file     = Channel.empty()
    nextclade_file  = mpx.out.nextclade_file
    vadr_file       = mpx.out.vadr_file
    dataset         = mpx.out.dataset
  } else if ( params.species == 'other') {
    other(fasta_prep.out.fastas.concat(multifastas).mix(cecret.out.consensus))
    pangolin_file   = Channel.empty()
    freyja_file     = Channel.empty()
    nextclade_file  = other.out.nextclade_file
    vadr_file       = other.out.vadr_file
    dataset         = other.out.dataset
  } else {
    pangolin_file   = Channel.empty()
    freyja_file     = Channel.empty()
    nextclade_file  = Channel.empty()
    vadr_file       = Channel.empty()
    dataset         = Channel.empty()
  }

  if ( params.relatedness ) { 
    msa(fasta_prep.out.fastas.concat(multifastas).concat(cecret.out.consensus), reference_genome, dataset) 

    tree      = msa.out.tree
    alignment = msa.out.msa
    matrix    = msa.out.matrix
  } else {
    tree      = Channel.empty()
    alignment = Channel.empty()
    matrix    = Channel.empty()
  }

  multiqc_combine(qc.out.fastqc_files.collect().ifEmpty([]),
    cecret.out.fastp_files.collect().ifEmpty([]),
    cecret.out.seqyclean_files1.collect().ifEmpty([]),
    cecret.out.seqyclean_files2.collect().ifEmpty([]),
    qc.out.kraken2_files.collect().ifEmpty([]),
    pangolin_file.collect().ifEmpty([]),
    cecret.out.ivar_files.collect().ifEmpty([]),
    qc.out.samtools_stats_files.collect().ifEmpty([]),
    qc.out.samtools_flagstat_files.collect().ifEmpty([]))

  cecret.out.consensus_results
    .mix(fasta_prep.out.fastas_results)
    // cecret subworkflow
    .join(cecret.out.cleaner_version,             remainder: true, by: 0 )
    .join(cecret.out.aligner_version,             remainder: true, by: 0 )
    .join(cecret.out.trimmer_version,             remainder: true, by: 0 )
    .join(cecret.out.ivar_version,                remainder: true, by: 0 )
    .join(cecret.out.fastp_results,               remainder: true, by: 0 )
    // qc subworkflow
    .join(qc.out.fastqc_1_results,                remainder: true, by: 0 )
    .join(qc.out.fastqc_2_results,                remainder: true, by: 0 )
    .join(qc.out.kraken2_target_results,          remainder: true, by: 0 )
    .join(qc.out.kraken2_human_results,           remainder: true, by: 0 )
    .join(qc.out.ivar_variants_results,           remainder: true, by: 0 )
    .join(qc.out.bcftools_variants_results,       remainder: true, by: 0 )
    .join(qc.out.insert_size_after_trimming,      remainder: true, by: 0 )
    .join(qc.out.samtools_coverage_results,       remainder: true, by: 0 )
    .join(qc.out.samtools_covdepth_results,       remainder: true, by: 0 )
    .join(qc.out.samtools_depth_results,          remainder: true, by: 0 )
    .join(qc.out.samtools_ampliconstats_results,  remainder: true, by: 0 )
    .join(qc.out.bedtools_results,                remainder: true, by: 0 )

    // seqyclean and anything from the organism-specific subworkflows will be added by pandas
    .set { results }

  summary(results)

  cecret.out.seqyclean_files1
    .collectFile(name: "Combined_SummaryStatistics.tsv",
      keepHeader: true,
      storeDir: "${params.outdir}/seqyclean")
    .set { seqyclean_file1 }

  cecret.out.seqyclean_files2
    .collectFile(name: "Combined_seqyclean_SummaryStatistics.tsv",
      keepHeader: true,
      storeDir: "${params.outdir}/seqyclean")
    .set { seqyclean_file2 }

  combine_results(nextclade_file.ifEmpty([]),
    pangolin_file.ifEmpty([]),
    vadr_file.ifEmpty([]),
    freyja_file.ifEmpty([]),
    seqyclean_file1.ifEmpty([]),
    seqyclean_file2.ifEmpty([]),
    summary.out.summary_file.collect().ifEmpty([]),
    combine_results_script)

  emit:
  bam       = cecret.out.bam_bai
  consensus = fasta_prep.out.fastas.mix(multifastas).mix(cecret.out.consensus)
  tree      = tree
  alignment = alignment
  matrix    = matrix
}

workflow.onComplete {
  println("Pipeline completed at: $workflow.complete")
  println("A summary of results can be found in a comma-delimited file: ${params.outdir}/summary/combined_summary.csv")
  println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
}
