process pango_collapse {
  tag        "SARS-CoV-2 lineage mapping"
  label      "process_low"
  publishDir "${params.outdir}", mode: 'copy'
  container  'quay.io/uphl/pango-collapse:0.7.2-2024-02-21'

  //#UPHLICA maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.pango_collapse

  input:
  file(file)

  output:
  path "pango_collapse/pango_collapse.csv", emit: results
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p pango_collapse logs/!{task.process}
    log=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log

    date > $log
    pango-collapse --version >> $log

    pango-collapse \
        !{params.pango_collapse_options} \
        --output pango_collapse/pango_collapse.csv \
        --lineage-column lineage \
        --collapse-file /pango-collapse/collapse.txt \
        !{file} \
        | tee -a $log
  '''
}
