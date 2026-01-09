process SNPDISTS {
  tag           "creating snp matrix with snp-dists"
  label         "process_medium"
  container     'staphb/snp-dists:1.2.0'

  input:
  file(msa)

  output:
  path "snp-dists/*.txt", emit: matrix
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: "${params.snpdists_options}"
  def prefix = task.ext.prefix ?: "snp-dists"
  """
    mkdir -p snp-dists logs/${task.process}
    log=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    date > \$log
    snp-dists -v >> \$log

    snp-dists ${args} \
      -j ${task.cpus} \
      ${msa} \
      > snp-dists/${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      snp-dists: \$(snp-dists -v | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}
