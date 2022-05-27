process fastqc {
  tag "${sample}"

  when:
  params.fastqc && sample != null

  input:
  tuple val(sample), file(fastq), val(type)

  output:
  path "fastqc/*.{html,zip}",                                             emit: files
  path "fastqc/*_fastqc.zip",                                             emit: fastqc_files
  tuple val(sample), env(raw_1),                                          emit: fastqc_1_results
  tuple val(sample), env(raw_2),                                          emit: fastqc_2_results
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}",  emit: log

  shell:
  '''
    mkdir -p fastqc logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    fastqc --version >> $log_file

    fastqc !{params.fastqc_options} \
      --outdir fastqc \
      --threads !{task.cpus} \
      !{fastq} \
      2>> $err_file >> $log_file

    zipped_fastq=($(ls fastqc/*fastqc.zip) "")

    raw_1=$(unzip -p ${zipped_fastq[0]} */fastqc_data.txt | grep "Total Sequences" | awk '{ print $3 }' )
    raw_2=NA
    if [ -f "${zipped_fastq[1]}" ] ; then raw_2=$(unzip -p fastqc/*fastqc.zip */fastqc_data.txt | grep "Total Sequences" | awk '{ print $3 }' ) ; fi

    if [ -z "$raw_1" ] ; then raw_1="0" ; fi
    if [ -z "$raw_2" ] ; then raw_2="0" ; fi
  '''
}
