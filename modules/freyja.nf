process freyja {
  tag "${sample}"
  label "medcpus"
  errorStrategy 'ignore'

  when:
  params.freyja

  input:
  tuple val(sample), file(bam), file(reference_genome)

  output:
  path "freyja/${sample}_demix.tsv",                                      emit: freyja_demix
  path "freyja/${sample}*",                                               emit: files
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{err,log}",  emit: log

  shell:
  '''
    mkdir -p freyja logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    # no version for 1.3.4

    freyja variants !{params.freyja_variants_options} \
      !{bam} \
      --variants freyja/!{sample}_variants.tsv \
      --depths freyja/!{sample}_depths.tsv \
      --ref !{reference_genome} \
      2>> $err_file >> $log_file

    freyja demix \
      !{params.freyja_demix_options} \
      freyja/!{sample}_variants.tsv \
      freyja/!{sample}_depths.tsv \
      --output freyja/!{sample}_demix.tsv \
      2>> $err_file >> $log_file

    freyja boot \
      freyja/!{sample}_variants.tsv \
      freyja/!{sample}_depths.tsv \
      --nt !{task.cpus} \
      --output_base freyja/!{sample}_boot.tsv \
      2>> $err_file >> $log_file
  '''
}

process freyja_aggregate {
  tag "Aggregating results from freyja"

  when:
  params.freyja_aggregate

  input:
  file(demix)

  output:
  path "freyja/aggregated*",                                                   emit: files
  path "freyja/aggregated-freyja.tsv",                                         emit: aggregated_freyja_file
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.{err,log}", emit: log

  shell:
  '''
    mkdir -p freyja logs/!{task.process} tmp
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null

    mv !{demix} tmp/.

    freyja aggregate !{params.freyja_aggregate_options} \
      tmp/ \
      --output freyja/aggregated-freyja.tsv \
      2>> $err_file >> $log_file

    freyja plot !{params.freyja_plot_options} \
      freyja/aggregated-freyja.tsv \
      --output freyja/aggregated-freyja.!{params.freyja_plot_filetype} \
      2>> $err_file >> $log_file
  '''
}
