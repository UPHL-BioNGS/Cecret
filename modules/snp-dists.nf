process snpdists {
  tag "createing snp matrix with snp-dists"

  when:
  params.snpdists

  input:
  file(msa)

  output:
  path "snp-dists/snp-dists.txt", emit: files
  path "logs/${task.process}/snp-dists.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p snp-dists logs/!{task.process}
    log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    date > $log
    snp-dists -v >> $log

    snp-dists !{params.snpdists_options} !{msa} > snp-dists/snp-dists.txt
  '''
}
