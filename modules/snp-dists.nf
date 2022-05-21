process snpdists {
  tag "createing snp matrix with snp-dists"

  when:
  params.snpdists

  input:
  file(msa)

  output:
  path "snp-dists/snp-dists.txt",                                         emit: files
  path "logs/${task.process}/snp-dists.${workflow.sessionId}.{log,err}",  emit: log

  shell:
  '''
    mkdir -p snp-dists logs/!{task.process}
    log_file=logs/!{task.process}/snp-dists.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/snp-dists.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    snp-dists -v >> $log_file

    snp-dists !{params.snpdists_options} !{msa} > snp-dists/snp-dists.txt 2> $err_file
  '''
}
