process nextclade {
  tag        "Clade Determination"
  label      "medcpus"
  publishDir "${params.outdir}", mode: 'copy'
  container  'nextstrain/nextclade:2.13.1'

  //#UPHLICA maxForks 10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.nextclade || params.msa == 'nextalign'

  input:
  file(fasta)

  output:
  path "nextclade/nextclade.csv",                                              emit: nextclade_file
  path "nextclade/*",                                                          emit: results
  path "nextclade/nextclade.aligned.fasta",                                    emit: nextclade_aligned_fasta
  path "dataset",                                                              emit: prepped_nextalign
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p nextclade dataset logs/!{task.process}
    log=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log

    date > $log
    nextclade --version >> $log
    nextclade_version=$(nextclade --version)

    nextclade dataset get --name !{params.nextclade_dataset} --output-dir dataset

    for fasta in !{fasta}
    do
      cat $fasta >> ultimate_fasta.fasta
    done

    nextclade run !{params.nextclade_options} \
      --input-dataset dataset \
      --output-all=nextclade/ \
      --jobs !{task.cpus} \
      ultimate_fasta.fasta \
      | tee -a $log
    cp ultimate_fasta.fasta nextclade/combined.fasta
  '''
}
