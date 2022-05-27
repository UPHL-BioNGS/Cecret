process pangolin {
  tag "SARS-CoV-2 lineage Determination"
  label "medcpus"

  when:
  params.pangolin

  input:
  file(fasta)

  output:
  path "pangolin/*",                                                            emit: results
  path "pangolin/lineage_report.csv",                                           emit: pangolin_file
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.{log,err}",  emit: log

  shell:
  '''
    mkdir -p pangolin logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    pangolin --all-versions >> $log_file

    for fasta in !{fasta}
    do
      cat $fasta >> ultimate_fasta.fasta
    done

    pangolin !{params.pangolin_options} \
      --threads !{task.cpus} \
      --outdir pangolin \
      ultimate_fasta.fasta \
      2>> $err_file >> $log_file
    cp ultimate_fasta.fasta pangolin/combined.fasta
  '''
}
