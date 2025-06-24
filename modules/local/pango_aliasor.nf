process PANGO_ALIASOR {
  tag        "SARS-CoV-2 lineage mapping"
  label      "process_low"
  container  'staphb/pango_aliasor:0.3.0-250625'

  input:
  file(file)

  output:
  path "pango_aliasor/*.csv", emit: results
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: "${params.pango_aliasor_options}"
  def prefix = task.ext.prefix ?: "pango_aliasor"
  """
    mkdir -p pango_aliasor logs/${task.process}
    log=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    date > \$log
    pip show pango_aliasor --version >> \$log

    aliasor.py \
      ${args} \
      --input ${file} \
      --output pango_aliasor/${prefix}.csv \
      --alias-key /key/alias_key.json \
      | tee -a \$log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      pango_aliasor: \$(pip show pango_aliasor | grep "Version" | awk '{print \$2}')
      container: ${task.container}
    END_VERSIONS
  """
}
