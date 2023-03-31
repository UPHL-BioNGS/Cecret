process fastqc {
  tag "${sample}"

  when:
  params.fastqc && sample != null

  input:
  tuple val(sample), file(fastq), val(type)

  output:
  path "fastqc/*.{html,zip}",                                             emit: files
  path "fastqc/*_fastqc.zip",                                             emit: fastqc_files
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p fastqc logs/!{task.process}
    log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log
    fastqc --version >> $log

    fastqc !{params.fastqc_options} \
      --outdir fastqc \
      --threads !{task.cpus} \
      !{fastq} \
      | tee -a $log
  '''
}
