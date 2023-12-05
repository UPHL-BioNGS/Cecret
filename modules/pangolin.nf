process pangolin {
  tag        "SARS-CoV-2 lineage Determination"
  label      "process_medium"
  publishDir "${params.outdir}", mode: 'copy'
  container  'staphb/pangolin:4.3.1-pdata-1.23.1-1'

  //#UPHLICA maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.pangolin

  input:
  file(fasta)

  output:
  path "pangolin/*",                                                            emit: results
  path "pangolin/lineage_report.csv",                                           emit: pangolin_file
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p pangolin logs/!{task.process}
    log=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log

    date > $log
    pangolin --all-versions >> $log

    for fasta in !{fasta}
    do
      cat $fasta >> ultimate_fasta.fasta
    done

    pangolin !{params.pangolin_options} \
      --threads !{task.cpus} \
      --outdir pangolin \
      --verbose \
      ultimate_fasta.fasta \
      | tee -a $log

    cp ultimate_fasta.fasta pangolin/combined.fasta
  '''
}
