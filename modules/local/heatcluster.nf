process HEATCLUSTER {
  tag        "HeatCluster"
  label      "process_low"
  container  'staphb/heatcluster:1.0.2c'

  input:
  file(matrix)

  output:
  path "heatcluster/*", optional : true, emit: for_multiqc
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: "${params.heatcluster_options}"
  def prefix = task.ext.prefix ?: "heatcluster"
  """
    mkdir -p heatcluster logs/${task.process}
    log_file=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log_file
    heatcluster.py -v >> \$log_file
    echo "container : ${task.container}" >> \$log_file
    echo "Nextflow command : " >> \$log_file
    cat .command.sh >> \$log_file

    heatcluster.py ${args} \
      -i ${matrix} \
      -o heatcluster/${prefix} \
      | tee -a \$log_file

    if [ -f "sorted_matrix.csv" ] ; then cp sorted_matrix.csv heatcluster/sorted_matrix.csv ; fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      heatcluster: \$(heatcluster.py -v | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}
