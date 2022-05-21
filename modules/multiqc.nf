process multiqc {
  tag "multiqc"

  when:
  params.multiqc

  input:
  file(fastqc)
  file(fastp)
  file(seqyclean)
  file(kraken2)
  file(pangolin)
  file(ivar)
  file(samtools_stats)
  file(samtools_flagstat)

  output:
  path "multiqc/multiqc_report.html",  optional: true,                         emit: html
  path "multiqc/multiqc_data/*",       optional: true,                         emit: files
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.{log,err}", emit: log

  shell:
  '''
    mkdir -p multiqc logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    multiqc --version >> $log_file

    multiqc !{params.multiqc_options} \
      --outdir multiqc \
      . \
      2>> $err_file >> $log_file
  '''
}
