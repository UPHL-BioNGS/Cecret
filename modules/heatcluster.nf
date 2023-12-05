process heatcluster {
  tag           "HeatCluster"
  publishDir    params.outdir, mode: 'copy'
  container     'quay.io/uphl/heatcluster:0.4.13-2023-12-05'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '10m'

  when:
  params.heatcluster

  input:
  file(matrix)

  output:
  path "heatcluster/heatcluster*"       , optional : true
  path "heatcluster/heatcluster_mqc.png", optional : true              , emit: for_multiqc
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log", emit: log_files

  shell:
  '''
    mkdir -p heatcluster logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    heatcluster.py -v >> $log_file
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    heatcluster.py !{params.heatcluster_options} \
        -i !{matrix} \
        -o heatcluster/heatcluster \
        | tee -a $log_file

    if [ -f "heatcluster/heatcluster.png" ] ; then cp heatcluster/heatcluster.png heatcluster/heatcluster_mqc.png ; fi
  '''
}
