process IQTREE {
  tag        "Creating phylogenetic tree with iqtree"
  label      "process_high"
  container  'staphb/iqtree3:3.0.1'

  input:
  file(msa)

  output:
  path "iqtree/*.{iqtree,treefile*,mldist,log}", emit: tree
  path "iqtree/*.treefile.nwk", emit: newick, optional: true
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: "${params.iqtree_options}"
  def prefix = task.ext.prefix ?: "iqtree"
  """
    mkdir -p iqtree logs/${task.process}
    log=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    date > \$log
    iqtree3 --version >> \$log

    if [ -n "${params.iqtree_outgroup}" ] && [ "${params.iqtree_outgroup}" != "null" ] && [ "${params.msa}" != "nextclade" ]
    then
      outgroup="-o ${params.iqtree_outgroup}"
      cat ${msa} | sed 's/${params.iqtree_outgroup}.*/${params.iqtree_outgroup}/g' > ${msa}.renamed
    else
      outgroup=""
      mv ${msa} ${msa}.renamed
    fi

    # creating a tree
    iqtree3 ${args} \
      -nt AUTO \
      -ntmax ${task.cpus} \
      -s ${msa}.renamed \
      -pre iqtree/${prefix} \
      \$outgroup \
      | tee -a \$log

    if [ -f "iqtree/iqtree.treefile" ] ; then cp iqtree/iqtree.treefile iqtree/iqtree.treefile.nwk ; fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      iqtree3: \$(iqtree3 --version | head -n 1 | awk '{print \$4}')
      container: ${task.container}
    END_VERSIONS
  """
}
