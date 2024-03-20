process heatcluster {
  tag           "HeatCluster"
  publishDir    path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/heatcluster:1.0.2c'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '10m'

  when:
  params.heatcluster && (task.ext.when == null || task.ext.when)

  input:
  file(matrix)

  output:
  path "heatcluster/*", optional : true, emit: for_multiqc
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log", emit: log_files
  path "versions.yml", emit: versions

  shell:
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
