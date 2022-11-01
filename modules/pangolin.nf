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
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p pangolin logs/!{task.process}
    log=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log

    date > $log
    pangolin --all-versions >> $log

    for fasta in !{fasta}
    do
      cat $fasta >> ultimate_fasta.fasta
    done

    pangolin !{params.pangolin_options} \
      --threads !{task.cpus} \
      --outdir pangolin \
      ultimate_fasta.fasta \
      | tee -a $log
    cp ultimate_fasta.fasta pangolin/combined.fasta
  '''
}
