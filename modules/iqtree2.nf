process iqtree2 {
  tag "Creating phylogenetic tree with iqtree"
  label "maxcpus"

  when:
  params.iqtree2

  input:
  file(msa)

  output:
  path "${task.process}/${task.process}.{iqtree,treefile,mldist,log}",         emit: files
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.{log,err}", emit: log

  shell:
  '''
    mkdir -p iqtree2 logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    iqtree2 --version >> $log_file

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
      >> $log_file 2>> $err_file
  '''
}
