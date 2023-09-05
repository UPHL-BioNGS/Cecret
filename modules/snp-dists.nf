process snpdists {
  tag        "creating snp matrix with snp-dists"
  label      "process_single"
  publishDir "${params.outdir}", mode: 'copy'
  container  'staphb/snp-dists:0.8.2'

  //#UPHLICA maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '45m'

  when:
  params.snpdists

  input:
  file(msa)

  output:
  path "snp-dists/snp-dists.txt", emit: matrix
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p snp-dists logs/!{task.process}
    log=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log

    date > $log
    snp-dists -v >> $log

    snp-dists !{params.snpdists_options} !{msa} > snp-dists/snp-dists.txt
  '''
}
