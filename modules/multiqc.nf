process multiqc_combine {
  tag "multiqc"

  when:
  params.multiqc

  input:
  file(fastqc)
  file(fastp)
  file(seqyclean)
  file(seqyclean2)
  file(kraken2)
  file(pangolin)
  file(ivar)
  file(samtools_stats)
  file(samtools_flagstat)

  output:
  path "multiqc/multiqc_report.html",  optional: true,                         emit: html
  path "multiqc/multiqc_data/*",       optional: true,                         emit: files
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p multiqc logs/!{task.process}
    log=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log
    multiqc --version >> $log

    multiqc !{params.multiqc_options} \
      --outdir multiqc \
      . \
      | tee -a $log
  '''
}
