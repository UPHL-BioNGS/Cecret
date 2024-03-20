process phytreeviz {
  tag           "Tree visualization"
  label         "maxcpus"
  publishDir    path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/phytreeviz:0.2.0'
  maxForks      10
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}

  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA cpus 14
  //#UPHLICA memory 60.GB
  //#UPHLICA time '24h'
  
  when:
  params.phytreeviz && (task.ext.when == null || task.ext.when)

  input:
  file(newick)

  output:
  path "phytreeviz/tree.png", emit: for_multiqc
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log", emit: log
  path "versions.yml", emit: versions

  shell:
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
