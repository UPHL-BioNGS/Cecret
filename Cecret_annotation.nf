#!/usr/bin/env nextflow

println("Currently using the Cecret workflow for use with amplicon-based Illumina hybrid library prep on MiSeq\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: v.20210208")
println("")

//# nextflow run Cecret/Cecret_annotation.nf -c Cecret/config/singularity.config

params.fastas = workflow.launchDir + '/fastas'
params.outdir = workflow.launchDir + '/cecret'
params.reference_genome = workflow.projectDir + "/configs/MN908947.3.fasta"

params.nextclade = true
params.pangolin = true

// for optional route of tree generation and counting snps between samples
params.relatedness = false
params.snpdists = true
params.iqtree = true
params.max_ambiguous = '0.50'
params.outgroup = 'MN908947.3'
params.mode='GTR'

params.maxcpus = Runtime.runtime.availableProcessors()
println("The maximum number of CPUS used in this workflow is ${params.maxcpus}")
if ( params.maxcpus < 5 ) {
  params.medcpus = params.maxcpus
} else {
  params.medcpus = 5
}

// This is where the results will be
println("The files and directory for results is " + params.outdir)

Channel
  .fromPath(params.reference_genome, type:'file')
  .ifEmpty{
    println("No reference genome was selected. Set with 'params.reference_genome'")
  }
  .set { reference_genome }

Channel
  .fromPath(["${params.fastas}/*.fa", "${params.fastas}/*.fasta"], type:'file')
  .ifEmpty{
    println("No fasta files were found. Set with 'params.fastas'.")
    exit 1
  }
  .into { fastas_mafft ; fastas_nextclade ; fastas_pangolin }

process pangolin {
  publishDir "${params.outdir}", mode: 'copy'
  tag "Lineage assignment with pangolin"
  echo false
  cpus 1

  when:
  params.pangolin

  input:
  file(fasta) from fastas_pangolin.collect()

  output:
  file("pangolin/lineage_report.csv")
  file("logs/pangolin/${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p pangolin logs/pangolin

    log_file=logs/pangolin/!{workflow.sessionId}.log
    err_file=logs/pangolin/!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    pangolin --version >> $log_file

    cat !{fasta} > ultimate.fasta

    pangolin --outdir pangolin/ ultimate.fasta 2>> $err_file >> $log_file
  '''
}

process nextclade {
  publishDir "${params.outdir}", mode: 'copy'
  tag "Clade assignment with nextclade"
  echo false
  cpus params.medcpus

  when:
  params.nextclade

  input:
  file(fasta) from fastas_nextclade.collect()

  output:
  file("nextclade/nextclade_report.tsv")
  file("logs/nextclade/${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p nextclade logs/nextclade
    log_file=logs/nextclade/!{workflow.sessionId}.log
    err_file=logs/nextclade/!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    nextclade --version >> $log_file

    cat !{fasta} > ultimate.fasta

    nextclade --jobs !{task.cpus} --input-fasta ultimate.fasta --output-tsv nextclade/nextclade_report.tsv 2>> $err_file >> $log_file
  '''
}

process mafft {
  publishDir "${params.outdir}", mode: 'copy'
  tag "Multiple Sequence Alignment with mafft"
  echo false
  cpus params.maxcpus

  input:
  file(fasta) from fastas_mafft.collect()
  file(reference) from reference_genome

  output:
  file("mafft/mafft_aligned.fasta") into msa_file
  file("mafft/mafft_aligned.fasta") into msa_file2
  file("logs/mafft/${workflow.sessionId}.{log,err}")

  when:
  params.relatedness

  shell:
  '''
    mkdir -p mafft logs/mafft
    log_file=logs/mafft/!{workflow.sessionId}.log
    err_file=logs/mafft/!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    echo "mafft version:" >> $log_file
    mafft --version >> $log_file

    cat !{reference}

    echo ">!{params.outgroup}" > reference.fasta
    grep -v ">" !{reference} >> reference.fasta

    cat !{fasta} > ultimate.fasta
    mafft --auto \
      --thread !{task.cpus} \
      --maxambiguous !{params.max_ambiguous} \
      --addfragments ultimate.fasta \
      reference.fasta \
      > mafft/mafft_aligned.fasta \
      2>> $err_file
  '''
}

process snpdists {
  publishDir "${params.outdir}", mode: 'copy'
  tag "SNP matrix with snp-dists"
  echo false
  cpus params.medcpus

  input:
  file(msa) from msa_file

  output:
  file("snp-dists/snp-dists.txt")
  file("logs/snp-dists/${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p snp-dists logs/snp-dists
    log_file=logs/snp-dists/!{workflow.sessionId}.log
    err_file=logs/snp-dists/!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    snp-dists -v >> $log_file

    snp-dists !{msa} > snp-dists/snp-dists.txt 2> $err_file
  '''
}

process iqtree {
  publishDir "${params.outdir}", mode: 'copy'
  tag "Creating phylogenetic tree with iqtree"
  echo false
  cpus params.maxcpus

  input:
  file(msa) from msa_file2

  output:
  file("iqtree/iqtree.{iqtree,treefile,mldist,log}")
  file("logs/iqtree/${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p iqtree logs/iqtree
    log_file=logs/iqtree/!{workflow.sessionId}.log
    err_file=logs/iqtree/!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    iqtree --version >> $log_file

    cat !{msa} | sed 's/!{params.outgroup}.*/!{params.outgroup}/g' > !{msa}.tmp
    mv !{msa}.tmp !{msa}

    # creating a tree
  	iqtree -ninit 2 \
      -n 2 \
      -me 0.05 \
      -nt AUTO \
      -ntmax !{task.cpus} \
      -s !{msa} \
      -pre iqtree/iqtree \
      -m !{params.mode} \
      -o !{params.outgroup} \
      >> $log_file 2>> $err_file
  '''
}

workflow.onComplete {
    println("Pipeline completed at: $workflow.complete")
    println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
}
