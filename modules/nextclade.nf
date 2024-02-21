process nextclade_dataset {
  tag        "Downloading NextClade Dataset"
  label      "process_medium"
  publishDir "${params.outdir}", mode: 'copy'
  container  'nextstrain/nextclade:3.2.1'

  //#UPHLICA maxForks 10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.nextclade

  output:
  path "dataset", emit: dataset
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p nextclade dataset logs/!{task.process}
    log=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log

    date > $log
    nextclade --version >> $log
    nextclade_version=$(nextclade --version)

    echo "Getting nextclade dataset for !{params.nextclade_dataset}" | tee -a $log
    nextclade dataset list | tee -a $log

    nextclade dataset get --name !{params.nextclade_dataset} --output-dir dataset
  '''
}

process nextclade {
  tag        "Clade Determination"
  label      "process_medium"
  publishDir "${params.outdir}", mode: 'copy'
  container  'nextstrain/nextclade:3.2.1'

  //#UPHLICA maxForks 10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.nextclade

  input:
  file(fasta)
  path(dataset)

  output:
  path "nextclade/nextclade.csv",                                                    emit: nextclade_file
  path "nextclade/*",                                                                emit: results
  tuple file("nextclade/nextclade.aligned.fasta"), file("nextclade/nextclade.nwk"), emit: prealigned, optional: true
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p nextclade dataset logs/!{task.process}
    log=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log

    date > $log
    nextclade --version >> $log
    nextclade_version=$(nextclade --version)

    for fasta in !{fasta}
    do
      cat $fasta >> ultimate_fasta.fasta
    done

    nextclade run !{params.nextclade_options} \
      --input-dataset !{dataset} \
      --output-all=nextclade/ \
      --jobs !{task.cpus} \
      ultimate_fasta.fasta \
      | tee -a $log
    cp ultimate_fasta.fasta nextclade/combined.fasta
  '''
}

process nextalign {
  tag        "Multiple Sequence Alignment"
  label      "process_medium"
  publishDir "${params.outdir}", mode: 'copy'
  container  'nextstrain/nextclade:3.2.1'

  //#UPHLICA maxForks 10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.nextclade

  input:
  file(consensus)
  path(dataset)

  output:
  path "nextalign/nextalign.aligned.fasta",   emit: msa
  path "nextalign/{*.fasta,nextalign.*.csv}", emit: files
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p nextalign logs/!{task.process}
    log=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log
    date > $log
    echo "nextalign version:" >> $log
    nextalign --version >> $log

    for fasta in !{consensus}
    do
      cat $fasta >> nextalign/ultimate.fasta
    done

    nextalign run !{params.nextalign_options} \
      --input-ref=!{dataset}/reference.fasta \
      --genemap=!{dataset}/genemap.gff \
      --jobs !{task.cpus} \
      --output-all=nextalign/ \
      nextalign/ultimate.fasta \
      | tee -a $log
  '''
}
