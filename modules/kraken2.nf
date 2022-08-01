process kraken2 {
  tag "${sample}"
  label "maxcpus"

  when:
  params.kraken2

  input:
  tuple val(sample), file(clean), val(paired_single), path(kraken2_db)

  output:
  path "kraken2/${sample}_kraken2_report.txt",                            emit: kraken2_files
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}",  emit: log
  tuple val(sample), env(percentage_cov),                                 emit: kraken2_target_results
  tuple val(sample), env(percentage_human),                               emit: kraken2_human_results

  shell:
  if ( paired_single == 'single' ) {
  '''
    mkdir -p kraken2 logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    kraken2 --version >> $log_file

    kraken2 !{params.kraken2_options} \
      --classified-out cseqs#.fq \
      --threads !{task.cpus} \
      --db !{kraken2_db} \
      !{clean} \
      --report kraken2/!{sample}_kraken2_report.txt \
      2>> $err_file >> $log_file

    percentage_human=$(grep "Homo sapiens" kraken2/!{sample}_kraken2_report.txt | awk '{if ($4 == "S") print $1}' | head -n 1)
    percentage_cov=$(grep "!{params.kraken2_organism}" kraken2/!{sample}_kraken2_report.txt | awk '{if ($4 == "S") print $1}' | head -n 1)

    if [ -z "$percentage_human" ] ; then percentage_human="0" ; fi
    if [ -z "$percentage_cov" ] ; then percentage_cov="0" ; fi
  '''
  } else if (paired_single == 'paired') {
    '''
      mkdir -p kraken2 logs/!{task.process}
      log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
      err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

      date | tee -a $log_file $err_file > /dev/null
      kraken2 --version >> $log_file

      kraken2 !{params.kraken2_options} \
        --paired \
        --classified-out cseqs#.fq \
        --threads !{task.cpus} \
        --db !{kraken2_db} \
        !{clean} \
        --report kraken2/!{sample}_kraken2_report.txt \
        2>> $err_file >> $log_file

      percentage_human=$(grep "Homo sapiens" kraken2/!{sample}_kraken2_report.txt | awk '{print $1}')
      percentage_cov=$(grep "!{params.kraken2_organism}" kraken2/!{sample}_kraken2_report.txt | awk '{print $1}')

      if [ -z "$percentage_human" ] ; then percentage_human="0" ; fi
      if [ -z "$percentage_cov" ] ; then percentage_cov="0" ; fi
    '''
  }
}
