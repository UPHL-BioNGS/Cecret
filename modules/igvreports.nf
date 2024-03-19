process igv_reports {
  tag         "${sample}"
  label       "process_high"
  publishDir  path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container   'quay.io/biocontainers/igv-reports:1.12.0--pyh7cba7a3_0'
  
  //#UPHLICA maxForks 10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.igv_reports && (task.ext.when == null || task.ext.when)

  input:
  tuple val(sample), file(vcf), file(bam), file(bai), file(reference_genome)

  output:
  path "igv_reports/*"
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"
  path "versions.yml", emit: versions

  shell:
  def args   = task.ext.args   ?: "${params.igv_reports_options}"
  def prefix = task.ext.prefix ?: "${sample}"
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
