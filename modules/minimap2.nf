process minimap2 {
  tag         "${sample}"
  label       "process_high"
  publishDir  path: "${params.outdir}", mode: 'copy', pattern: 'logs/*/*log'
  container   'staphb/minimap2:2.26'

  //#UPHLICA maxForks 10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  task.ext.when == null || task.ext.when

  input:
  tuple val(sample), file(reads), file(reference_genome)

  output:
  tuple val(sample), file("aligned/${sample}.sam"), emit: sam
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"
  tuple val("${params.aligner}"), env(minimap2_version), emit: aligner_version
  path "versions.yml", emit: versions

  shell:
  def args   = task.ext.args   ?: "${params.minimap2_options}}"
  def fastq  = reads.join(" ")
  def prefix = task.ext.prefix ?: "${sample}"
  """
    mkdir -p aligned logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log
    minimap2 --version >> \$log
    minimap2_version=\$(echo "minimap2 : "\$(minimap2 --version))

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
