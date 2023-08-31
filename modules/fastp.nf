process fastp {
  tag        "${sample}"
  label      "process_single"
  publishDir "${params.outdir}", mode: 'copy'
  container  'staphb/fastp:0.23.4'

  //#UPHLICA maxForks 10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '45m'

  when:
  sample != null

  input:
  tuple val(sample), file(reads), val(paired_single)

  output:
  tuple val(sample), file("fastp/${sample}_{clean_PE1,clean_PE2,cln}.fastq.gz"), optional: true,  emit: clean_reads
  path "fastp/${sample}_fastp.html",                                                              emit: html
  tuple val(sample), file("fastp/${sample}_fastp.json"),                                          emit: fastp_files
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}"
  tuple val("${params.cleaner}"), env(cleaner_version),                                           emit: cleaner_version

  shell:
  if ( paired_single == 'paired' ) {
    '''
      mkdir -p fastp logs/!{task.process}
      log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
      err=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

      # time stamp + capturing tool versions
      date > $log
      fastp --version >> $log
      cleaner_version="$(fastp --version 2>&1 | head -n 1)"

      fastp !{params.fastp_options} \
        -i !{reads[0]} \
        -I !{reads[1]} \
        -o fastp/!{sample}_clean_PE1.fastq.gz \
        -O fastp/!{sample}_clean_PE2.fastq.gz \
        -h fastp/!{sample}_fastp.html \
        -j fastp/!{sample}_fastp.json \
        2>> $err | tee -a $log
    '''
  } else if ( paired_single == 'single' ) {
    '''
      mkdir -p fastp logs/!{task.process}
      log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
      err=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

      # time stamp + capturing tool versions
      date > $log
      fastp --version >> $log
      cleaner_version="$(fastp --version 2>&1 | head -n 1)"

      fastp !{params.fastp_options} \
        -i !{reads} \
        -o fastp/!{sample}_cln.fastq.gz \
        -h fastp/!{sample}_fastp.html \
        -j fastp/!{sample}_fastp.json \
        2>> $err | tee -a $log
    '''
  }
}
