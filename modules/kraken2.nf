process kraken2 {
  tag "${sample}"
  label "maxcpus"

  when:
  params.kraken2

  input:
  tuple val(sample), file(clean), path(kraken2_db)

  output:
  path "kraken2/${sample}_kraken2_report.txt",                            emit: kraken2_files
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"

  shell:
  if ( clean =~ "cln" ) {
  '''
    mkdir -p kraken2 logs/!{task.process}
    log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    date > $log
    kraken2 --version >> $log

    kraken2 !{params.kraken2_options} \
      --classified-out cseqs#.fq \
      --threads !{task.cpus} \
      --db !{kraken2_db} \
      !{clean} \
      --report kraken2/!{sample}_kraken2_report.txt \
      | tee -a $log

  '''
  } else {
    '''
      mkdir -p kraken2 logs/!{task.process}
      log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

      date > $log
      kraken2 --version >> $log

      kraken2 !{params.kraken2_options} \
        --paired \
        --classified-out cseqs#.fq \
        --threads !{task.cpus} \
        --db !{kraken2_db} \
        !{clean} \
        --report kraken2/!{sample}_kraken2_report.txt \
        | tee -a $log
    '''
  }
}
