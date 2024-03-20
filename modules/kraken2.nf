process kraken2 {
  tag        "${sample}"
  label      "process_high"
  publishDir path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container  'staphb/kraken2:2.1.3'

  //#UPHLICA maxForks 10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.kraken2 && (task.ext.when == null || task.ext.when)

  input:
  tuple val(sample), file(clean), path(kraken2_db)

  output:
  path "kraken2/${sample}_kraken2_report.txt", emit: kraken2_files
  path "kraken2/*", emit: everything
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"
  path "versions.yml", emit: versions

  shell:
  def args   = task.ext.args   ?: "${params.kraken2_options}"
  def prefix = task.ext.prefix ?: "${sample}"
  def reads  = clean.join(" ")
  if ( clean =~ "cln" ) {
  """
    mkdir -p kraken2 logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    date > \$log
    kraken2 --version >> \$log

    kraken2 ${args} \
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
    END_VERSIONS
  """
  } else {
  """
    mkdir -p kraken2 logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    date > \$log
    kraken2 --version >> \$log

    kraken2 ${args} \
      --paired \
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
      container: ${task.container}
    END_VERSIONS
  """
  }
}
