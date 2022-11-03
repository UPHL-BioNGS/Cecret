process bwa {
  tag "${sample}"
  label "maxcpus"

  input:
  tuple val(sample), file(reads), file(reference_genome)

  output:
  tuple val(sample), file("bwa/${sample}.sam"),                           emit: sam
  tuple val(sample), env(bwa_version),                                    emit: aligner_version
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p bwa logs/!{task.process}
    log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log
    echo "bwa $(bwa 2>&1 | grep Version )" >> $log
    bwa_version="bwa : "$(bwa 2>&1 | grep Version)

    # index the reference fasta file
    bwa index !{reference_genome}

    # bwa mem command
    bwa mem -t !{task.cpus} !{reference_genome} !{reads} > bwa/!{sample}.sam
  '''
}
