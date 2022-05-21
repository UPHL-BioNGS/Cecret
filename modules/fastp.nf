process fastp {
  tag "${sample}"

  when:
  sample != null

  input:
  tuple val(sample), file(reads), val(paired_single)

  output:
  tuple val(sample), file("fastp/${sample}_clean_PE{1,2}.fastq.gz"),              optional: true, emit: paired_files
  tuple val(sample), file("fastp/${sample}_cln.fastq.gz"),                        optional: true, emit: single_files
  typle val(sample), file("fastp/${sample}_{clean_PE1,clean_PE2,cln}.fastq.gz"),  optional: true, emit: clean_type
  path "fastp/${sample}_fastp.html",                                                              emit: html
  path "fastp/${sample}_fastp.json",                                                              emit: fastp_files
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}",                          emit: log
  tuple val(sample), env(passed_reads),                                                           emit: fastp_results
  tuple val(sample), env(cleaner_version),                                                        emit: cleaner_version

  shell:
  if ( paired_single == 'paired' ) {
    '''
      mkdir -p fastp logs/!{task.process}
      log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
      err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

      # time stamp + capturing tool versions
      date | tee -a $log_file $err_file > /dev/null
      fastp --version >> $log_file 2>> $err_file
      cleaner_version="$(fastp --version 2>&1 | head -n 1)"

      fastp !{params.fastp_options} \
        -i !{reads[0]} \
        -I !{reads[1]} \
        -o !{task.process}/!{sample}_clean_PE1.fastq.gz \
        -O !{task.process}/!{sample}_clean_PE2.fastq.gz \
        -h !{task.process}/!{sample}_fastp.html \
        -j !{task.process}/!{sample}_fastp.json \
        2>> $err_file >> $log_file

      passed_reads=$(grep "reads passed filter" $err_file | tail -n 1 | cut -f 2 -d ":" | sed 's/ //g' )
      if [ -z "$passed_reads" ] ; then passed_reads="0" ; fi
    '''
  } else if ( paired_single == 'single' ) {
    '''
      mkdir -p fastp logs/!{task.process}
      log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
      err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

      # time stamp + capturing tool versions
      date | tee -a $log_file $err_file > /dev/null
      fastp --version >> $log_file 2>> $err_file
      cleaner_version="$(fastp --version 2>&1 | head -n 1)"

      fastp !{params.fastp_options} \
        -i !{reads} \
        -o !{task.process}/!{sample}_cln.fastq.gz \
        -h !{task.process}/!{sample}_fastp.html \
        -j !{task.process}/!{sample}_fastp.json \
        2>> $err_file >> $log_file

      passed_reads=$(grep "reads passed filter" $err_file | tail -n 1 | cut -f 2 -d ":" | sed 's/ //g' )
      if [ -z "$passed_reads" ] ; then passed_reads="0" ; fi
    '''
  }
}
