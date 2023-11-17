process iqtree2 {
  tag        "Creating phylogenetic tree with iqtree"
  label      "process_high"
  publishDir "${params.outdir}", mode: 'copy'
  container  'staphb/iqtree2:2.2.2.6'

  //#UPHLICA maxForks 10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.iqtree2

  input:
  file(msa)

  output:
  path "iqtree2/iqtree2.{iqtree,treefile*,mldist,log}", emit: tree
  path "iqtree2/iqtree2.treefile.nwk"                 , emit: newick, optional: true
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p iqtree2 logs/!{task.process}
    log=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log

    date > $log
    iqtree2 --version >> $log

    if [ -n "!{params.iqtree2_outgroup}" ] && [ "!{params.iqtree2_outgroup}" != "null" ] && [ "!{params.msa}" != "nextclade" ]
    then
      outgroup="-o !{params.iqtree2_outgroup}"
      cat !{msa} | sed 's/!{params.iqtree2_outgroup}.*/!{params.iqtree2_outgroup}/g' > !{msa}.renamed
    else
      outgroup=""
      mv !{msa} !{msa}.renamed
    fi

    # creating a tree
    iqtree2 !{params.iqtree2_options} \
      -nt AUTO \
      -ntmax !{task.cpus} \
      -s !{msa}.renamed \
      -pre iqtree2/iqtree2 \
      $outgroup \
      | tee -a $log

    if [ -f "iqtree2/iqtree2.treefile" ] ; then cp iqtree2/iqtree2.treefile iqtree2/iqtree2.treefile.nwk ; fi
  '''
}
