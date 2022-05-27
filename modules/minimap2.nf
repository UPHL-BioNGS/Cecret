process minimap2 {
  tag "${sample}"
  label "maxcpus"

  input:
  tuple val(sample), file(reads), file(reference_genome)

  output:
  tuple val(sample), file("aligned/${sample}.sam"),                       emit: sam
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}",  emit: log
  tuple val(sample), env(minimap2_version),                               emit: aligner_version

  shell:
  '''
    mkdir -p aligned logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    minimap2 --version >> $log_file
    minimap2_version=$(echo "minimap2 : "$(minimap2 --version))

    minimap2 !{params.minimap2_options} \
      -ax sr -t !{task.cpus} \
      -o aligned/!{sample}.sam \
      !{reference_genome} !{reads} 2>> $err_file >> $log_file
  '''
}
