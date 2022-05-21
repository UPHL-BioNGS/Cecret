process mafft {
  tag "Multiple Sequence Alignment"
  label "maxcpus"
  errorStrategy 'retry'
  maxRetries 2

  input:
  tuple file(fasta), file(reference_genome)

  output:
  path "mafft/mafft_aligned.fasta",                                   emit: msa
  path "logs/${task.process}/mafft.${workflow.sessionId}.{log,err}",  emit: log

  shell:
  '''
    mkdir -p mafft logs/!{task.process}
    log_file=logs/!{task.process}/mafft.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/mafft.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    echo "mafft version:" >> $log_file
    mafft --version 2>&1 >> $log_file

    for fasta in !{fasta}
    do
      cat $fasta >> mafft/ultimate.fasta
    done

    mafft --auto \
      !{params.mafft_options} \
      --thread !{task.cpus} \
      --addfragments mafft/ultimate.fasta \
      !{reference_genome} \
      > mafft/mafft_aligned.fasta \
      2>> $err_file
  '''
}
