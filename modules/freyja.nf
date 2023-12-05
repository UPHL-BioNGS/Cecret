process freyja_variants {
  tag           "${sample}"
  label         "process_medium"
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  publishDir    "${params.outdir}", mode: 'copy'
  container     'quay.io/uphl/freyja:1.4.7-2023-12-05'

  //#UPHLICA maxForks 10
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.freyja

  input:
  tuple val(sample), file(bam), file(reference_genome)

  output:
  tuple val(sample), file("freyja/${sample}_{depths,variants}.tsv"), optional: true, emit: variants
  path "freyja/${sample}*",                                          optional: true, emit: files
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p freyja logs/!{task.process}
    log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    date > $log
    freyja --version >> $log

    freyja variants !{params.freyja_variants_options} \
      !{bam} \
      --variants freyja/!{sample}_variants.tsv \
      --depths freyja/!{sample}_depths.tsv \
      --ref !{reference_genome} \
      | tee -a $log
  '''
}

process freyja_demix {
  tag           "${sample}"
  label         "process_medium"
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  publishDir    "${params.outdir}", mode: 'copy'
  container     'quay.io/uphl/freyja:1.4.7-2023-12-05'

  //#UPHLICA maxForks 10
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.freyja

  input:
  tuple val(sample), file(variants)

  output:
  path "freyja/${sample}_demix.tsv", optional: true, emit: demix
  path "freyja/${sample}*",          optional: true, emit: files
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p freyja logs/!{task.process}
    log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    date > $log
    freyja --version >> $log

    freyja demix \
      !{params.freyja_demix_options} \
      !{variants[1]} \
      !{variants[0]} \
      --output freyja/!{sample}_demix.tsv \
      | tee -a $log
  '''
}

process freyja_boot {
  tag           "${sample}"
  label         "process_medium"
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  publishDir    "${params.outdir}", mode: 'copy'
  container     'quay.io/uphl/freyja:1.4.7-2023-12-05'

  //#UPHLICA maxForks 10
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.freyja

  input:
  tuple val(sample), file(variants)

  output:
  path "freyja/${sample}*", emit: files
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p freyja logs/!{task.process}
    log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    date > $log
    freyja --version >> $log

    freyja boot !{params.freyja_boot_options} \
      !{variants[1]} \
      !{variants[0]} \
      --nt !{task.cpus} \
      --output_base freyja/!{sample}_boot.tsv \
      | tee -a $log
  '''
}

process freyja_aggregate {
  tag        "Aggregating results from freyja"
  label      "process_single"
  publishDir "${params.outdir}", mode: 'copy'
  container  'quay.io/uphl/freyja:1.4.7-2023-12-05'

  //#UPHLICA maxForks 10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.freyja_aggregate

  input:
  file(demix)
  file(graphs)

  output:
  path "freyja/aggregated*",                                                   emit: files
  path "freyja/aggregated-freyja.tsv",                                         emit: aggregated_freyja_file
  path "freyja/*mqc.png",                                                      emit: for_multiqc
  path "freyja/*png"
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p freyja logs/!{task.process} tmp
    log=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log

    date > $log

    mv !{demix} tmp/.

    freyja aggregate !{params.freyja_aggregate_options} \
      tmp/ \
      --output freyja/aggregated-freyja.tsv \
      | tee -a $log

    freyja plot !{params.freyja_plot_options} \
      freyja/aggregated-freyja.tsv \
      --output freyja/aggregated-freyja.!{params.freyja_plot_filetype} \
      | tee -a $log

    python3 !{graphs}
  '''
}
