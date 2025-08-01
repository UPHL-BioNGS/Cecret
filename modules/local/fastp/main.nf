process FASTP {
  tag        "${meta.id}"
  label      "process_single"
  container  'staphb/fastp:1.0.1'

  input:
  tuple val(meta), file(reads)

  output:
  tuple val(meta), file("fastp/*_{clean_PE1,clean_PE2,cln}.fastq.gz"), optional: true,  emit: clean_reads
  path "fastp/*_fastp.html", emit: html
  tuple val(meta), file("fastp/*_fastp.json"), emit: fastp_files
  path "logs/${task.process}/*", emit: log
  tuple val("${params.cleaner}"), env("cleaner_version"), emit: cleaner_version
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: "${params.fastp_options}"
  def prefix = task.ext.prefix ?: "${meta.id}"
  if ( meta.single_end ) {
    """
      mkdir -p fastp logs/${task.process}
      log=logs/${task.process}/${prefix}.${workflow.sessionId}.log
      err=logs/${task.process}/${prefix}.${workflow.sessionId}.err

      # time stamp + capturing tool versions
      date > \$log
      fastp --version >> \$log
      cleaner_version=\$(fastp --version 2>&1 | awk '{print \$NF}')

      fastp ${args} \
        -i ${reads} \
        -o fastp/${prefix}_cln.fastq.gz \
        -h fastp/${prefix}_fastp.html \
        -j fastp/${prefix}_fastp.json \
        2>> \$err | tee -a \$log

      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
        fastp: \$(fastp --version 2>&1 | awk '{print \$NF}')
        container: ${task.container}
      END_VERSIONS
    """
  } else {
    """
    mkdir -p fastp logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log
    err=logs/${task.process}/${prefix}.${workflow.sessionId}.err

    # time stamp + capturing tool versions
    date > \$log
    fastp --version >> \$log
    cleaner_version=\$(fastp --version 2>&1 | awk '{print \$NF}')

    fastp ${args} \
      -i ${reads[0]} \
      -I ${reads[1]} \
      -o fastp/${prefix}_clean_PE1.fastq.gz \
      -O fastp/${prefix}_clean_PE2.fastq.gz \
      -h fastp/${prefix}_fastp.html \
      -j fastp/${prefix}_fastp.json \
      2>> \$err | tee -a \$log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      fastp: \$(fastp --version 2>&1 | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
    """
  }
}
