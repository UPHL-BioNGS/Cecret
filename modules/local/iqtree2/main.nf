process IQTREE2 {
  tag        "Creating phylogenetic tree with iqtree"
  label      "process_high"
  container  'staphb/iqtree2:2.4.0'

  input:
  file(msa)

  output:
  path "iqtree2/iqtree2.{iqtree,treefile*,mldist,log}", emit: tree
  path "iqtree2/iqtree2.treefile.nwk", emit: newick, optional: true
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: "${params.iqtree2_options}"
  def prefix = task.ext.prefix ?: "iqtree2"
  """
    mkdir -p iqtree2 logs/${task.process}
    log=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    date > \$log
    iqtree2 --version >> \$log

    if [ -n "${params.iqtree2_outgroup}" ] && [ "${params.iqtree2_outgroup}" != "null" ] && [ "${params.msa}" != "nextclade" ]
    then
      outgroup="-o ${params.iqtree2_outgroup}"
      cat ${msa} | sed 's/${params.iqtree2_outgroup}.*/${params.iqtree2_outgroup}/g' > ${msa}.renamed
    else
      outgroup=""
      mv ${msa} ${msa}.renamed
    fi

    # creating a tree
    iqtree2 ${args} \
      -nt AUTO \
      -ntmax ${task.cpus} \
      -s ${msa}.renamed \
      -pre iqtree2/${prefix} \
      \$outgroup \
      | tee -a \$log

    if [ -f "iqtree2/iqtree2.treefile" ] ; then cp iqtree2/iqtree2.treefile iqtree2/iqtree2.treefile.nwk ; fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      iqtree2: \$(iqtree2 --version | head -n 1 | awk '{print \$4}')
      container: ${task.container}
    END_VERSIONS
  """
}
