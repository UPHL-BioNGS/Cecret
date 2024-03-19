process fastqc {
  tag        "${sample}"
  label      "process_single"
  publishDir params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container  'staphb/fastqc:0.12.1'

  //#UPHLICA maxForks 10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '45m'

  when:
  params.fastqc && sample != null && (task.ext.when == null || task.ext.when)

  input:
  tuple val(sample), file(fastq), val(type)

  output:
  path "fastqc/*.html",                   emit: files
  path "fastqc/*_fastqc.zip",             emit: fastqc_files
  path "fastqc/${sample}_fastq_name.csv", emit: fastq_name
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"
  path "versions.yml", emit: versions

  shell:
  def args   = task.ext.args   ?: "${params.fastqc_options}"
  def reads  = fastq.join(" ")
  def prefix = task.ext.prefix ?: "${sample}"
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
