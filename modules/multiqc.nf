process multiqc_combine {
  tag        "multiqc"
  label      "process_single"
  publishDir "${params.outdir}", mode: 'copy'
  container  'quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0'

  //#UPHLICA maxForks 10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '45m'

  when:
  params.multiqc

  input:
  file(input)

  output:
  path "multiqc/multiqc_report.html",  optional: true,                    emit: html
  path "multiqc/multiqc_data/*",       optional: true,                    emit: files
  path "multiqc/multiqc_data",         optional: true,                    emit: multiqc_data
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
