process PANGOLIN {
  tag        "SARS-CoV-2 lineage Determination"
  label      "process_medium"
  container  'staphb/pangolin:4.3.4-pdata-1.37'

  input:
  file(fasta)

  output:
  path "pangolin/*", emit: results
  path "pangolin/lineage_report.csv", emit: pangolin_file
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
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
