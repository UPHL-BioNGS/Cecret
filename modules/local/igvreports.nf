process IGV_REPORTS {
  tag         "${meta.id}"
  label       "process_high"
  container   'staphb/igv-reports:1.12.0'

  input:
  tuple val(meta), file(vcf), file(bam), file(bai), file(reference_genome)

  output:
  val(meta), emit: meta // for linter
  path "igv_reports/*", emit: report
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: "${params.igv_reports_options}"
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p igv_reports logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log
    echo "igv-reports does not print version to screen" >> \$log

    create_report ${args} \
      --fasta ${reference_genome} \
      --tracks ${vcf} ${bam} \
      --output igv_reports/${prefix}_igvjs_viewer.html \
      ${vcf} \
      | tee -a \$log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      igv_reports: NT
      container: ${task.container}
    END_VERSIONS
  """
}
