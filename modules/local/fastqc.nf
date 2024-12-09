process FASTQC {
  tag        "${meta.id}"
  label      "process_single"
  container  'staphb/fastqc:0.12.1'

  input:
  tuple val(meta), file(fastq)

  output:
  val(meta), emit: meta // for linter
  path "fastqc/*.html", emit: files
  path "fastqc/*_fastqc.zip", emit: fastqc_files
  path "fastqc/*_fastq_name.csv", emit: fastq_name
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: "${params.fastqc_options}"
  def reads  = fastq.join(" ")
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p fastqc logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log
    fastqc --version >> \$log

    fastqc ${args} \
      --outdir fastqc \
      --threads ${task.cpus} \
      ${reads} \
      | tee -a \$log

    echo "sample,fastq_1,fastq_2"             > fastqc/${prefix}_fastq_name.csv
    echo "${prefix},${fastq[0]},${fastq[1]}" >> fastqc/${prefix}_fastq_name.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      fastqc: \$(fastqc --version | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}
