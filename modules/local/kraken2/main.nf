process KRAKEN2 {
  tag        "${meta.id}"
  label      "process_high"
  container  'staphb/kraken2:2.1.6-viral-20250402'

  input:
  tuple val(meta), file(clean)

  output:
  val(meta), emit: meta // for linter
  path "kraken2/*_kraken2_report.txt", emit: kraken2_files
  path "kraken2/*", emit: everything
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: "${params.kraken2_options}"
  def prefix = task.ext.prefix ?: "${meta.id}"
  def reads  = clean.join(" ")
  def paired = meta.single_end ? "" : "--paired"
  """
    mkdir -p kraken2 logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    date > \$log
    kraken2 --version >> \$log

    kraken2 ${args} \
      ${paired} \
      --classified-out kraken2/${prefix}.cseqs#.fastq.gz \
      --unclassified-out kraken2/${prefix}.useqs#.fastq.gz \
      --threads ${task.cpus} \
      --db /kraken2-db \
      ${reads} \
      --report kraken2/${prefix}_kraken2_report.txt \
      | tee -a \$log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      kraken2: \$(kraken2 --version | head -n 1 | awk '{print \$NF}')
    END_VERSIONS
  """
}

process KRAKEN2_DB {
  tag        "${meta.id}"
  label      "process_high"
  container  'staphb/kraken2:2.1.6'

  input:
  tuple val(meta), file(clean), path(kraken2_db)

  output:
  val(meta), emit: meta // for linter
  path "kraken2/*_kraken2_report.txt", emit: kraken2_files
  path "kraken2/*", emit: everything
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: "${params.kraken2_options}"
  def prefix = task.ext.prefix ?: "${meta.id}"
  def reads  = clean.join(" ")
  def paired = meta.single_end ? "" : "--paired"
  """
    mkdir -p kraken2 logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    date > \$log
    kraken2 --version >> \$log

    kraken2 ${args} \
      ${paired} \
      --classified-out kraken2/${prefix}.cseqs#.fastq.gz \
      --unclassified-out kraken2/${prefix}.useqs#.fastq.gz \
      --threads ${task.cpus} \
      --db ${kraken2_db} \
      ${reads} \
      --report kraken2/${prefix}_kraken2_report.txt \
      | tee -a \$log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      kraken2: \$(kraken2 --version | head -n 1 | awk '{print \$NF}')
      kraken2_db : ${params.kraken2_db}
    END_VERSIONS
  """
}