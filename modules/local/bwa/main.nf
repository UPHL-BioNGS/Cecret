process BWA {
  tag         "${meta.id}"
  label       "process_high"
  container   'staphb/bwa:0.7.19'

  input:
  tuple val(meta), file(reads), file(reference_genome)

  output:
  tuple val(meta), file("bwa/*.sam"), emit: sam
  tuple val("${params.aligner}"), env("bwa_version"), emit: aligner_version
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: ''
  def fastq  = reads.join(" ")
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p bwa logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log
    echo "bwa \$(bwa 2>&1 | grep Version )" >> \$log
    bwa_version=\$(bwa 2>&1 | grep Version | awk '{print \$NF}')

    # index the reference fasta file
    bwa index ${reference_genome}

    # bwa mem command
    bwa mem ${args} \
      -t ${task.cpus} \
      ${reference_genome} \
      ${fastq} \
      > bwa/${prefix}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      bwa: \$(bwa 2>&1 | grep Version | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}
