process phytreeviz {
  tag           "Tree visualization"
  label         "maxcpus"
  publishDir    params.outdir, mode: 'copy'
  container     'quay.io/uphl/phytreeviz:0.1.0-2023-11-15'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA cpus 14
  //#UPHLICA memory 60.GB
  //#UPHLICA time '24h'
  
  input:
  file(newick)

  output:
  path "phytreeviz/tree.png"
  path "phytreeviz/tree_mqc.png"                                       , emit: for_multiqc
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log", emit: log

  shell:
  '''
    mkdir -p phytreeviz logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    phytreeviz -v >> $log_file
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    phytreeviz !{params.phytreeviz_options} \
        -i !{newick} \
        -o phytreeviz/tree.png

    if [ -f "phytreeviz/tree.png" ] ; then cp phytreeviz/tree.png phytreeviz/tree_mqc.png ; fi 
  '''
}
