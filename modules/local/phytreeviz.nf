process PHYTREEVIZ {
  tag           "Tree visualization"
  label         "process_medium"
  container     'staphb/phytreeviz:0.2.0'
  
  
  input:
  file(newick)

  output:
  path "phytreeviz/tree.png", emit: for_multiqc
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: "${params.phytreeviz_options}"
  def prefix = task.ext.prefix ?: "tree"
  """
    mkdir -p phytreeviz logs/${task.process}
    log_file=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log_file
    phytreeviz -v >> \$log_file
    echo "container : ${task.container}" >> \$log_file
    echo "Nextflow command : " >> \$log_file
    cat .command.sh >> \$log_file

    phytreeviz ${args} \
      -i ${newick} \
      -o phytreeviz/${prefix}.png \
      | tee -a \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      phytreeviz: \$(phytreeviz -v | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}
