#!/usr/bin/env nextflow

println("For annotating SARS-CoV-2 fastas with pangolin, nextclade, and vadr\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: v.0.20210611")
println("")

params.fastas = workflow.launchDir + '/fastas'
params.outdir = workflow.launchDir + '/cecret'
params.reference_genome = workflow.projectDir + "/configs/MN908947.3.fasta"

params.maxcpus = Runtime.runtime.availableProcessors()
println("The maximum number of CPUS used in this workflow is ${params.maxcpus}")
if ( params.maxcpus < 5 ) {
  params.medcpus = params.maxcpus
} else {
  params.medcpus = 5
}

println("The files and directory for results is " + params.outdir)

Channel
  .fromPath(["${params.fastas}/*.fa", "${params.fastas}/*.fasta"], type:'file')
  .ifEmpty{
    println("No fasta files were found. Set directory with 'params.fastas'.")
    exit 1
  }
  .set { fastas }

process fasta_prep {
  publishDir "${params.outdir}", mode: 'copy', overwrite: true
  tag "${fasta}"
  echo false
  cpus 1
  container 'staphb/parallel-perl:latest'

  input:
  file(fasta) from fastas

  output:
  file("${task.process}/${fasta}") into fastas_mafft, fastas_pangolin, fastas_nextclade, fastas_vadr

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{workflow.sessionId}.err

    name=$(echo !{fasta} | sed 's/.fa.*//g')

    echo ">$name" > !{task.process}/!{fasta}
    grep -v ">" !{fasta} | fold -w 75 >> !{task.process}/!{fasta}
  '''
}

params.pangolin = true
params.pangolin_options = ''
process pangolin {
  publishDir "${params.outdir}", mode: 'copy'
  tag "Lineage assignment with pangolin"
  echo false
  cpus 1
  container 'staphb/pangolin:latest'

  when:
  params.pangolin

  input:
  file(fasta) from fastas_pangolin.collect()

  output:
  file("${task.process}/lineage_report.csv")
  file("${task.process}/ultimate.fasta")
  file("logs/${task.process}/${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    pangolin --version >> $log_file

    cat !{fasta} > !{task.process}/ultimate.fasta

    pangolin !{params.pangolin_options} \
      --outdir !{task.process} \
      !{task.process}/ultimate.fasta \
      2>> $err_file >> $log_file
  '''
}

params.nextclade = true
params.nextclade_options = ''
process nextclade {
  publishDir "${params.outdir}", mode: 'copy'
  tag "Clade assignment with nextclade"
  echo false
  cpus params.medcpus
  container 'neherlab/nextclade:latest'

  when:
  params.nextclade

  input:
  file(fasta) from fastas_nextclade.collect()

  output:
  file("${task.process}/nextclade_report.tsv")
  file("${task.process}/ultimate.fasta")
  file("logs/${task.process}/${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    nextclade --version >> $log_file

    cat !{fasta} > !{task.process}/ultimate.fasta

    nextclade !{params.nextclade_options} \
      --jobs !{task.cpus} \
      --input-fasta !{task.process}/ultimate.fasta \
      --output-tsv !{task.process}/nextclade_report.tsv \
      2>> $err_file >> $log_file
  '''
}

params.vadr = true
if ( Math.round(Runtime.runtime.totalMemory() / 10241024) / 2 > params.medcpus && params.vadr ) {
  vadrmemory = params.medcpus + params.medcpus
  vadrcpus = params.medcpus
} else {
  vadrmemory = 2
  vadrcpus = 1
}

params.vadr_options = '--split --glsearch -s  -r --nomisc --lowsim5term 2 --lowsim3term 2 --alt_fail lowscore,fstukcnf,insertnn,deletinn'
params.vadr_reference = 'sarscov2'
params.vadr_mdir = '/opt/vadr/vadr-models'
process vadr {
  publishDir "${params.outdir}", mode: 'copy'
  tag "Fasta QC with vadr"
  echo false
  cpus vadrcpus
  memory vadrmemory.GB
  container 'staphb/vadr:latest'
  stageInMode = 'symlink'

  when:
  params.vadr

  input:
  file(fasta) from fastas_vadr.collect()

  output:
  file("${task.process}/*")
  file("logs/${task.process}/${task.process}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    echo "no version" >> $log_file
    v-annotate.pl -h >> $log_file

    cat !{fasta} > ultimate.fasta

    v-annotate.pl !{params.vadr_options} \
      --noseqnamemax \
      --mkey !{params.vadr_reference} \
      --mdir !{params.vadr_mdir} \
      --cpu !{task.cpus} \
      ultimate.fasta \
      !{task.process} \
      2>> $err_file >> $log_file
  '''
}

params.relatedness = false
if (params.relatedness){
  Channel
    .fromPath(params.reference_genome, type:'file')
    .ifEmpty{
      println("No reference genome was selected. Set with 'params.reference_genome'")
    }
    .set { reference_genome }

  params.mafft_options = ''
  process mafft {
    publishDir "${params.outdir}", mode: 'copy'
    tag "Multiple Sequence Alignment with mafft"
    echo false
    cpus params.maxcpus
    container 'staphb/mafft:latest'
    errorStrategy 'retry'
    maxRetries 3

    input:
    file(fasta) from fastas_mafft.collect()
    file(reference) from reference_genome

    output:
    file("${task.process}/mafft_aligned.fasta") into msa_file
    file("${task.process}/mafft_aligned.fasta") into msa_file2
    file("logs/${task.process}/${workflow.sessionId}.{log,err}")

    when:
    params.relatedness

    shell:
    '''
      mkdir -p !{task.process} logs/!{task.process}
      log_file=logs/!{task.process}/!{workflow.sessionId}.log
      err_file=logs/!{task.process}/!{workflow.sessionId}.err

      date | tee -a $log_file $err_file > /dev/null
      echo "mafft version:" >> $log_file
      mafft --version >> $log_file

      echo ">!{params.outgroup}" > reference.fasta
      grep -v ">" !{reference} >> reference.fasta

      cat !{fasta} > !{task.process}/ultimate.fasta
      mafft !{params.mafft_options} \
        --auto \
        --thread !{task.cpus} \
        --maxambiguous !{params.max_ambiguous} \
        --addfragments !{task.process}/ultimate.fasta \
        reference.fasta \
        > !{task.process}/mafft_aligned.fasta \
        2>> $err_file
    '''
  }

  params.snpdists = true
  params.snpdists_options = ''
  process snpdists {
    publishDir "${params.outdir}", mode: 'copy'
    tag "SNP matrix with snp-dists"
    echo false
    cpus params.medcpus
    container 'staphb/snp-dists:latest'

    input:
    file(msa) from msa_file

    when:
    params.snpdists

    output:
    file("snp-dists/snp-dists.txt")
    file("logs/${task.process}/${workflow.sessionId}.{log,err}")

    shell:
    '''
      mkdir -p snp-dists logs/!{task.process}
      log_file=logs/!{task.process}/!{workflow.sessionId}.log
      err_file=logs/!{task.process}/!{workflow.sessionId}.err

      date | tee -a $log_file $err_file > /dev/null
      snp-dists -v >> $log_file

      snp-dists !{params.snpdists_options} !{msa} > snp-dists/snp-dists.txt 2> $err_file
    '''
  }

  params.iqtree_options = ''
  params.iqtree = true
  params.max_ambiguous = '0.50'
  params.outgroup = 'MN908947.3'
  params.mode='GTR'
  process iqtree {
    publishDir "${params.outdir}", mode: 'copy'
    tag "Creating phylogenetic tree with iqtree"
    echo false
    cpus params.maxcpus
    container 'staphb/iqtree:latest'

    input:
    file(msa) from msa_file2

    when:
    params.iqtree

    output:
    file("${task.process}/iqtree.{iqtree,treefile,mldist,log}")
    file("logs/${task.process}/${workflow.sessionId}.{log,err}")

    shell:
    '''
      mkdir -p !{task.process} logs/!{task.process}
      log_file=logs/!{task.process}/!{workflow.sessionId}.log
      err_file=logs/!{task.process}/!{workflow.sessionId}.err

      date | tee -a $log_file $err_file > /dev/null
      iqtree --version >> $log_file

      cat !{msa} | sed 's/!{params.outgroup}.*/!{params.outgroup}/g' > !{msa}.tmp
      mv !{msa}.tmp !{msa}

      # creating a tree
    	iqtree !{params.iqtree_options} \
        -ninit 2 \
        -n 2 \
        -me 0.05 \
        -nt AUTO \
        -ntmax !{task.cpus} \
        -s !{msa} \
        -pre !{task.process}/iqtree \
        -m !{params.mode} \
        -o !{params.outgroup} \
        >> $log_file 2>> $err_file
    '''
  }
}

workflow.onComplete {
    println("Pipeline completed at: $workflow.complete")
    println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
}
