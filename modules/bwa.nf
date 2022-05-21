process bwa {
  tag "${sample}"
  label "maxcpus"

  input:
  tuple val(sample), file(reads), file(reference_genome)

  output:
  tuple val(sample), file("bwa/${sample}.sam"),                           emit: sam
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}",  emit: log
  tuple val(sample), env(bwa_version),                                    emit: aligner_version

  shell:
  '''
    mkdir -p bwa logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "bwa $(bwa 2>&1 | grep Version )" >> $log_file
    bwa_version="bwa : "$(bwa 2>&1 | grep Version)

    # index the reference fasta file
    bwa index !{reference_genome}

    # bwa mem command
    bwa mem -t !{task.cpus} !{reference_genome} !{reads} 2>> $err_file > bwa/!{sample}.sam
  '''
}
