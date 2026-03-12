process HEATCLUSTER {
  tag           "HeatCluster"
  label         "process_single"
  container     'staphb/heatcluster:1.3.0'

  input:
  file(matrix)

  output:
  path "heatcluster/*", optional : true, emit: files
  path "heatcluster/*.png", optional : true, emit: for_multiqc
  path "logs/${task.process}/*.log", emit: log_files
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: '--cluster-t 10'
  def prefix = task.ext.prefix ?: "heatcluster"
  """
    mkdir -p heatcluster tmp logs/${task.process}
    log_file=logs/${task.process}/heatcluster.${workflow.sessionId}.log

    export MPLCONFIGDIR=tmp

    heatcluster ${args} \
      -i ${matrix} \
      -o heatcluster/${prefix}.png  \
      -c heatcluster/${prefix}_sorted.csv \
      | tee -a \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      heatcluster: \$(echo \$(heatcluster -v | awk '{print \$NF}' ))
      container: ${task.container}
    END_VERSIONS
  """
}