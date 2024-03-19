process iqtree2 {
  tag        "Creating phylogenetic tree with iqtree"
  label      "process_high"
  publishDir path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container  'staphb/iqtree2:2.2.2.7'

  //#UPHLICA maxForks 10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.iqtree2 && (task.ext.when == null || task.ext.when)

  input:
  file(msa)

  output:
  path "iqtree2/iqtree2.{iqtree,treefile*,mldist,log}", emit: tree
  path "iqtree2/iqtree2.treefile.nwk", emit: newick, optional: true
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log"
  path "versions.yml", emit: versions

  shell:
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
