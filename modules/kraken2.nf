process kraken2 {
  tag        "${sample}"
  label      "process_high"
  publishDir "${params.outdir}", mode: 'copy'
  container  'staphb/kraken2:2.1.3'

  //#UPHLICA maxForks 10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.kraken2

  input:
  tuple val(sample), file(clean), path(kraken2_db)

  output:
  path "kraken2/${sample}_kraken2_report.txt",                            emit: kraken2_files
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"

  shell:
  if ( clean =~ "cln" ) {
  '''
    mkdir -p kraken2 logs/!{task.process}
    log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    date > $log
    kraken2 --version >> $log

    kraken2 !{params.kraken2_options} \
      --classified-out cseqs#.fq \
      --threads !{task.cpus} \
      --db !{kraken2_db} \
      !{clean} \
      --report kraken2/!{sample}_kraken2_report.txt \
      | tee -a $log

  '''
  } else {
    '''
      mkdir -p kraken2 logs/!{task.process}
      log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

      date > $log
      kraken2 --version >> $log

      kraken2 !{params.kraken2_options} \
        --paired \
        --classified-out cseqs#.fq \
        --threads !{task.cpus} \
        --db !{kraken2_db} \
        !{clean} \
        --report kraken2/!{sample}_kraken2_report.txt \
        | tee -a $log
    '''
  }
}
