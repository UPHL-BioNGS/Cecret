#!/usr/bin/env nextflow

println("Currently using the Cecret workflow for use with amplicon-based Illumina hybrid library prep on MiSeq\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: v.3.0.20220520")
println("")

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

//# default reference files for SARS-CoV-2 (part of the github repository)
params.reference_genome                     = workflow.projectDir + "/configs/MN908947.3.fasta"
params.gff_file                             = workflow.projectDir + "/configs/MN908947.3.gff"
params.primer_set                           = 'ncov_V4'
if ( params.primer_set == 'ncov_V3' ) {
  params.primer_bed                         = workflow.projectDir + "/configs/artic_V3_nCoV-2019.primer.bed"
  params.amplicon_bed                       = workflow.projectDir + "/configs/artic_V3_nCoV-2019.insert.bed"
} else if ( params.primer_set == 'ncov_V4' ) {
  params.primer_bed                         = workflow.projectDir + "/configs/artic_V4_SARS-CoV-2.primer.bed"
  params.amplicon_bed                       = workflow.projectDir + "/configs/artic_V4_SARS-CoV-2.insert.bed"
} else if ( params.primer_set == 'ncov_V4.1' ) {
  params.primer_bed                         = workflow.projectDir + "/configs/artic_V4.1_SARS-CoV-2.primer.bed"
  params.amplicon_bed                       = workflow.projectDir + "/configs/artic_V4.1_SARS-CoV-2.insert.bed"
} else {
  println("!{params.primer_set} has not been defined as an acceptable value for 'params.primer_set'.")
  println("Current acceptable values are" )
  println("SARS-CoV-2 artic primer V3 : 'params.primer_set' = 'ncov_V3'" )
  println("SARS-CoV-2 artic primer V4 : 'params.primer_set' = 'ncov_V4'" )
  println("SARS-CoV-2 artic primer V4.1 (Version 4 with spike in) : 'params.primer_set' = 'ncov_V4.1'" )
  exit 1
}

//# specifying the core workflow
params.trimmer                              = 'ivar'
params.cleaner                              = 'seqyclean'
params.aligner                              = 'bwa'
params.msa                                  = 'mafft'

//# to toggle off processes
params.bcftools_variants                    = true // YOLO
params.fastqc                               = true
params.ivar_variants                        = true
params.samtools_stats                       = true
params.samtools_coverage                    = true
params.samtools_depth                       = true
params.samtools_flagstat                    = true
params.samtools_ampliconstats               = true
params.samtools_plot_ampliconstats          = true
params.bedtools_multicov                    = true
params.nextclade                            = true
params.pangolin                             = true
params.kraken2                              = false
params.filter                               = false
params.vadr                                 = true
params.freyja                               = true
params.freyja_aggregate                     = true
params.multiqc                              = true

//# for optional route of tree generation and counting snps between samples
params.relatedness                          = false
params.snpdists                             = true
params.iqtree2                              = true

//# parameters for processes
params.fastqc_options                       = ''
params.seqyclean_contaminant_file           = "/Adapters_plus_PhiX_174.fasta"
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
params.pangolin_options                     = ''
params.nextclade_options                    = ''
params.nextclade_dataset                    = 'sars-cov-2'
params.vadr_options                         = '--split --glsearch -s -r --nomisc --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn'
params.vadr_reference                       = 'sarscov2'
params.vadr_mdir                            = '/opt/vadr/vadr-models'
params.vadr_trim_options                    = '--minlen 50 --maxlen 30000'
//# for optional contamination determination
params.kraken2_db                           = false
params.kraken2_organism                     = "Severe acute respiratory syndrome-related coronavirus"
params.freyja_variants_options              = ''
params.freyja_demix_options                 = ''
params.freyja_boot_options                  = '--nb 1000'
params.freyja_aggregate_options             = ''
params.freyja_plot_options                  = ''
params.freyja_plot_filetype                 = 'png'
params.mafft_options                        = '--maxambiguous 0.5'
params.nextalign_options                    = '--genes E,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S --include-reference'
params.snpdists_options                     = '-csv'
params.iqtree2_outgroup                     = 'MN908947'
params.iqtree2_options                      = '-ninit 2 -n 2 -me 0.05 -m GTR'
params.multiqc_options                      = ''

include { fasta_prep ; summary; combine_results } from './modules/cecret.nf'      addParams(fastqc: params.fastqc,
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
                                                                                            ivar_consensus_options: params.ivar_consensus_options)
include { qc }                                    from './subworkflows/qc'        addParams(fastqc: params.fastqc,
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
                                                                                            samtools_plot_ampliconstats_options: params.samtools_plot_ampliconstats_options,
                                                                                            freyja: params.freyja,
                                                                                            freyja_variants_options: params.freyja_variants_options,
                                                                                            freyja_demix_options: params.freyja_demix_options,
                                                                                            freyja_aggregate: params.freyja_aggregate,
                                                                                            freyja_aggregate_options: params.freyja_aggregate_options,
                                                                                            freyja_plot_options: params.freyja_plot_options)
include { annotation }                            from './subworkflows/annotation' addParams(vadr: params.vadr,
                                                                                            vadr_options: params.vadr_options,
                                                                                            vadr_reference: params.vadr_reference,
                                                                                            vadr_mdir: params.vadr_mdir,
                                                                                            pangolin: params.pangolin,
                                                                                            pangolin_options: params.pangolin_options,
                                                                                            nextclade: params.nextclade,
                                                                                            nextclade_options: params.nextclade_options,
                                                                                            nextclade_dataset: params.nextclade_dataset)
include { msa }                                   from './subworkflows/msa'       addParams(msa: params.msa,
                                                                                            nextalign_options: params.nextalign_options,
                                                                                            mafft_options: params.mafft_options,
                                                                                            iqtree2: params.iqtree2,
                                                                                            iqtree2_options: params.iqtree2_options,
                                                                                            iqtree2_outgroup: params.iqtree2_outgroup,
                                                                                            snpdists: params.snpdists,
                                                                                            snpdists_options: params.snpdists_options)
include { multiqc }                               from './modules/multiqc'        addParams(multiqc: params.multiqc,
                                                                                            multiqc_options: params.multiqc_options)

//# getting input files
Channel
  .fromFilePairs(["${params.reads}/*_R{1,2}*.{fastq,fastq.gz,fq,fq.gz}",
                  "${params.reads}/*{1,2}*.{fastq,fastq.gz,fq,fq.gz}"], size: 2 )
  .map { reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1], "paired" ) }
  .view { "Paired-end Fastq files found : ${it[0]}" }
  .set { paired_reads }

Channel
  .fromPath("${params.single_reads}/*.{fastq,fastq.gz,fq,fq.gz}")
  .map { reads -> tuple(reads.simpleName, reads, "single" ) }
  .view { "Fastq files found : ${it[0]}" }
  .set { single_reads }

Channel
  .fromPath("${params.fastas}/*{.fa,.fasta,.fna}", type:'file')
  .map { fasta -> tuple(fasta.baseName, fasta ) }
  .view { "Fasta file found : ${it[0]}" }
  .set { fastas }

Channel
  .fromPath("${params.multifastas}/*{.fa,.fasta,.fna}", type:'file')
  .view { "MultiFasta file found : ${it}" }
  .set { multifastas }

//# Checking for input files and giving an explanatory message if none are found
paired_reads
  .concat(single_reads)
  .concat(fastas)
  .concat(multifastas)
  .ifEmpty{
    println("FATAL : No input files were found!")
    println("No paired-end fastq files were found at ${params.reads}. Set 'params.reads' to directory with paired-end reads")
    println("No single-end fastq files were found at ${params.single_reads}. Set 'params.single_reads' to directory with single-end reads")
    println("No fasta files were found at ${params.fastas}. Set 'params.fastas' to directory with fastas.")
    println("No multifasta files were found at ${params.fastas}. Set 'params.multifastas' to directory with multifastas.")
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

Channel
  .fromPath(params.gff_file, type:'file')
  .view { "GFF file for Reference Genome : $it"}
  .set { gff_file }

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

kraken2_db = params.kraken2_db
  ? Channel.fromPath(params.kraken2_db, type:'dir').view { "Kraken2 database : $it" }
  : Channel.empty()

//# getting scripts
Channel
  .fromPath("${workflow.projectDir}/bin/combine_results.py", type:'file')
  .set { combine_results_script }

// This is where the results will be
println("The files and directory for results is " + params.outdir)
println("A table summarizing results will be created: ${params.outdir}/cecret_results.csv\n")

paired_reads
  .concat(single_reads)
  .ifEmpty{ println("No fastq or fastq.gz files were found at ${params.reads} or ${params.single_reads}") }
  .set { reads }

workflow {
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

    annotation(fasta_prep.out.fastas.concat(multifastas).concat(cecret.out.consensus))

    if ( params.relatedness ) { msa(fasta_prep.out.fastas.concat(multifastas).concat(cecret.out.consensus), reference_genome, annotation.out.dataset) }

    //multiqc(qc.out.fastqc_files, cecret.out.fastp_files)
    // , cecret.out.seqyclean_files.collect(),
    //   qc.out.kraken2_files.collect(),
    //   annotation.out.pangolin.collect(),
    //   cecret.out.ivar_files.collect(),
    //   qc.out.samtools_stats_files.collect(),
    //   qc.out.samtools_flagstat_files.collect())

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
      .join(qc.out.samtools_covdepth_results,       remainder: true, by: 0)
      .join(qc.out.samtools_depth_results,          remainder: true, by: 0 )
      .join(qc.out.samtools_ampliconstats_results,  remainder: true, by: 0 )
      .join(qc.out.bedtools_results,                remainder: true, by: 0 )
      // seqyclean and anything from the annotation subworkflow will be added by pandas
      .set { results }

      summary(results)

      cecret.out.seqyclean_files1
        .collectFile(name: "Combined_SummaryStatistics.tsv",
          keepHeader: true,
          sort: true,
          storeDir: "${params.outdir}/seqyclean")
        .set { seqyclean_file1 }

      cecret.out.seqyclean_files2
        .collectFile(name: "Combined_seqyclean_SummaryStatistics.tsv",
          keepHeader: true,
          sort: true,
          storeDir: "${params.outdir}/seqyclean")
        .set { seqyclean_file2 }

      combine_results(annotation.out.nextclade_file.ifEmpty([]),
        annotation.out.pangolin_file.ifEmpty([]),
        annotation.out.vadr_file.ifEmpty([]),
        qc.out.freyja_file.ifEmpty([]),
        seqyclean_file1.ifEmpty([]),
        seqyclean_file2.ifEmpty([]),
        summary.out.summary_file.collect().ifEmpty([]),
        combine_results_script)
}

workflow.onComplete {
  println("Pipeline completed at: $workflow.complete")
  println("A summary of results can be found in a comma-delimited file: ${params.outdir}/summary/combined_summary.csv")
  println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
}
