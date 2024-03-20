process pangolin {
  tag        "SARS-CoV-2 lineage Determination"
  label      "process_medium"
  publishDir path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container  'staphb/pangolin:4.3.1-pdata-1.26'

  //#UPHLICA maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.pangolin && (task.ext.when == null || task.ext.when)

  input:
  file(fasta)

  output:
  path "pangolin/*", emit: results
  path "pangolin/lineage_report.csv", emit: pangolin_file
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log"
  path "versions.yml", emit: versions

  shell:
  def args   = task.ext.args   ?: "${params.pangolin_options}"
  def file   = fasta.join(" ")
  def prefix = task.ext.prefix ?: "combined"
  """
    mkdir -p pangolin logs/${task.process}
    log=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    date > \$log
    pangolin --all-versions >> \$log

    for fasta in ${file}
    do
      cat \$fasta >> ultimate_fasta.fasta
    done

    pangolin ${args} \
      --threads ${task.cpus} \
      --outdir pangolin \
      --verbose \
      ultimate_fasta.fasta \
      | tee -a \$log

    cp ultimate_fasta.fasta pangolin/${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    \$(pangolin --all-versions | awk '{(\$1 = \$1 ); print "  " \$1 ": " \$2}' | sed 's/::/:/g')
      container: ${task.container}
    END_VERSIONS
  """
}
