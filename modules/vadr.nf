process vadr {
  tag "QC metrics"
  label "medcpus"

  when:
  params.vadr

  input:
  file(fasta)

  output:
  path "${task.process}/*",             optional: true,                         emit: vadr_files
  path "${task.process}/vadr.vadr.sqa", optional: true,                         emit: vadr_file
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.{log,err}",  emit: log

  shell:
  '''
    mkdir -p logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    echo "no version" >> $log_file
    v-annotate.pl -h >> $log_file

    for fasta in !{fasta}
    do
      cat $fasta >> ultimate_fasta.fasta
    done

    fasta-trim-terminal-ambigs.pl !{params.vadr_trim_options} \
      ultimate_fasta.fasta > trimmed_ultimate_fasta.fasta

    if [ -s "trimmed_ultimate_fasta.fasta" ]
    then
      v-annotate.pl !{params.vadr_options} \
        --cpu !{task.cpus} \
        --noseqnamemax \
        --mkey !{params.vadr_reference} \
        --mdir !{params.vadr_mdir} \
        trimmed_ultimate_fasta.fasta \
        vadr \
        2>> $err_file >> $log_file
    fi
    cp ultimate_fasta.fasta vadr/combined.fasta
    cp trimmed_ultimate_fasta.fasta vadr/trimmed.fasta
  '''
}
