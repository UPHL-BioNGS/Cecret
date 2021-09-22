#!/usr/bin/env nextflow

println("For annotating SARS-CoV-2 fastas with pangolin, nextclade, and vadr\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: v.1.2.20210930")
println("")

params.fastas = workflow.launchDir + '/fastas'
params.outdir = workflow.launchDir + '/cecret_annotation'
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

params.vadr = true
params.pangolin = true
params.relatedness = false
params.nextclade = true
if (params.nextclade) {
  Channel
    .fromPath(params.reference_genome, type:'file')
    .ifEmpty{
      println("No reference genome was selected. Set with 'params.reference_genome'")
    }
    .set { reference_genome_nextclade }
}

process fasta_prep {
  publishDir "${params.outdir}", mode: 'copy', overwrite: true
  tag "${fasta}"
  cpus 1
  container 'staphb/parallel-perl:latest'

  input:
  file(fasta) from fastas

  output:
  file("${task.process}/${fasta}") into prepped_fastas, fastas_mafft, fastas_pangolin, fastas_nextclade, fastas_vadr

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

prepped_fastas
  .collectFile(name: "Ultimate.fasta", storeDir: "${params.outdir}")
  .into { multifasta_pangolin ; multifasta_vadr ; multifasta_nextclade ; multifasta_mafft }

params.pangolin_options = ''
process pangolin {
  publishDir "${params.outdir}", mode: 'copy'
  tag "Lineage assignment with pangolin"
  cpus 1
  container 'staphb/pangolin:latest'

  when:
  params.pangolin

  input:
  file(fasta) from multifasta_pangolin

  output:
  file("${task.process}/lineage_report.csv")
  file("logs/${task.process}/${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    pangolin --version >> $log_file

    pangolin !{params.pangolin_options} \
      --outdir !{task.process} \
      !{fasta} \
      2>> $err_file >> $log_file
  '''
}

params.nextclade_options = ''
params.nextclade_genes = 'E,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S'
process nextclade {
  publishDir "${params.outdir}", mode: 'copy'
  tag "Clade assignment with nextclade"
  cpus params.medcpus
  container 'nextstrain/nextclade:latest'

  when:
  params.nextclade

  input:
  file(fasta) from multifasta_nextclade
  file(reference) from reference_genome_nextclade

  output:
  file("${task.process}/*")
  file("logs/${task.process}/${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    nextclade --version >> $log_file

    wget https://raw.githubusercontent.com/nextstrain/nextclade/master/data/sars-cov-2/genemap.gff
    wget https://raw.githubusercontent.com/nextstrain/nextclade/master/data/sars-cov-2/tree.json
    wget https://raw.githubusercontent.com/nextstrain/nextclade/master/data/sars-cov-2/qc.json
    wget https://raw.githubusercontent.com/nextstrain/nextclade/master/data/sars-cov-2/primers.csv

    nextclade !{params.nextclade_options} \
      --input-fasta=!{fasta} \
      --input-root-seq=!{reference} \
      --genes=!{params.nextclade_genes} \
      --input-gene-map=genemap.gff \
      --input-tree=tree.json \
      --input-qc-config=qc.json \
      --input-pcr-primers=primers.csv \
      --output-json=!{task.process}/nextclade.json \
      --output-csv=!{task.process}/nextclade.csv \
      --output-tsv=!{task.process}/nextclade.tsv \
      --output-tree=!{task.process}/nextclade.auspice.json \
      --output-dir=!{task.process} \
      --output-basename=!{task.process} \
      2>> $err_file >> $log_file
  '''
}

params.vadr_options = '--split --glsearch -s -r --nomisc --mkey sarscov2 --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn'
params.vadr_reference = 'sarscov2'
params.vadr_mdir = '/opt/vadr/vadr-models'
params.vadr_trim_options = '--minlen 50 --maxlen 30000'
process vadr {
  publishDir "${params.outdir}", mode: 'copy'
  tag "Fasta QC with vadr"
  cpus params.medcpus
  container 'staphb/vadr:latest'

  when:
  params.vadr

  input:
  file(fasta) from multifasta_vadr

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

    fasta-trim-terminal-ambigs.pl !{params.vadr_trim_options} \
      !{fasta} > trimmed_!{fasta}

    v-annotate.pl !{params.vadr_options} \
      --cpu !{task.cpus} \
      --noseqnamemax \
      --mkey !{params.vadr_reference} \
      --mdir !{params.vadr_mdir} \
      trimmed_!{fasta} \
      !{task.process} \
      2>> $err_file >> $log_file
  '''
}

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
    cpus params.maxcpus
    container 'staphb/mafft:latest'
    errorStrategy 'retry'
    maxRetries 3

    input:
    file(fasta) from multifasta_mafft
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

      mafft !{params.mafft_options} \
        --auto \
        --thread !{task.cpus} \
        --maxambiguous !{params.max_ambiguous} \
        --addfragments !{fasta} \
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

  params.iqtree2_options = ''
  params.iqtree2 = true
  params.max_ambiguous = '0.50'
  params.outgroup = 'MN908947.3'
  params.mode='GTR'
  process iqtree2 {
    publishDir "${params.outdir}", mode: 'copy'
    tag "Creating phylogenetic tree with iqtree2"
    cpus params.maxcpus
    container 'staphb/iqtree2:latest'

    input:
    file(msa) from msa_file2

    when:
    params.iqtree2

    output:
    file("${task.process}/iqtree2.{iqtree,treefile,mldist,log}")
    file("logs/${task.process}/${workflow.sessionId}.{log,err}")

    shell:
    '''
      mkdir -p !{task.process} logs/!{task.process}
      log_file=logs/!{task.process}/!{workflow.sessionId}.log
      err_file=logs/!{task.process}/!{workflow.sessionId}.err

      date | tee -a $log_file $err_file > /dev/null
      iqtree2 --version >> $log_file

      cat !{msa} | sed 's/!{params.outgroup}.*/!{params.outgroup}/g' > !{msa}.tmp
      mv !{msa}.tmp !{msa}

      # creating a tree
    	iqtree2 !{params.iqtree2_options} \
        -ninit 2 \
        -n 2 \
        -me 0.05 \
        -nt AUTO \
        -ntmax !{task.cpus} \
        -s !{msa} \
        -pre !{task.process}/iqtree2 \
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
