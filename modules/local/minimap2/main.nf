process MINIMAP2 {
  tag         "${meta.id}"
  label       "process_high"
  container   'staphb/minimap2:2.30'

  input:
  tuple val(meta), file(reads), file(reference_genome)

  output:
  tuple val(meta), file("aligned/*.sam"), emit: sam
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: "${params.minimap2_options}"
  def fastq  = reads.join(" ")
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p aligned logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log
    minimap2 --version >> \$log
    minimap2_version=\$(minimap2 --version | awk '{print \$NF}')

    minimap2 ${args} \
      -ax sr -t ${task.cpus} \
      -o aligned/${prefix}.sam \
      ${reference_genome} ${fastq} \
      | tee -a \$log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      minimap2: \$(minimap2 --version | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}
