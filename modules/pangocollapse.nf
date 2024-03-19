process pango_collapse {
  tag        "SARS-CoV-2 lineage mapping"
  label      "process_low"
  publishDir path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container  'quay.io/uphl/pango-collapse:0.8.2-2024-03-19'

  //#UPHLICA maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.pango_collapse && (task.ext.when == null || task.ext.when)

  input:
  file(file)

  output:
  path "pango_collapse/pango_collapse.csv", emit: results
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log"
  path "versions.yml", emit: versions

  shell:
  def args   = task.ext.args   ?: "${params.pango_collapse_options}"
  def prefix = task.ext.prefix ?: "pango_collapse"
  """
    mkdir -p pango_collapse logs/${task.process}
    log=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    date > \$log
    pango-collapse --version >> \$log

    pango-collapse \
        ${args} \
        --output pango_collapse/${prefix}.csv \
        --lineage-column lineage \
        --collapse-file /pango-collapse/collapse.txt \
        ${file} \
        | tee -a \$log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      pango-collapse: \$(pango-collapse --version | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}
