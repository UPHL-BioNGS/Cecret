process seqyclean {
  tag "${sample}"

  when:
  sample != null

  input:
  tuple val(sample), file(reads), val(paired_single)

  output:
  tuple val(sample), file("seqyclean/${sample}_clean_PE{1,2}.fastq.gz"),                                    optional: true, emit: paired_reads
  tuple val(sample), file("seqyclean/${sample}_cln_SE.fastq.gz"),                                           optional: true, emit: single_reads
  tuple val(sample), file("seqyclean/${sample}_{cln_SE,clean_PE1,clean_PE2}.fastq.gz"), val(paired_single), optional: true, emit: clean_reads
  path "seqyclean/${sample}_cl*n_SummaryStatistics.tsv",                                                                    emit: seqyclean_files_collect
  path "seqyclean/${sample}_cl*n_SummaryStatistics.txt",                                                                    emit: txt
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}",                                                    emit: log
  tuple val(sample), env(cleaner_version),                                                                                  emit: cleaner_version

  shell:
  if ( paired_single == "single" ) {
  '''
    mkdir -p seqyclean logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "seqyclean version: $(seqyclean -h | grep Version)" >> $log_file
    cleaner_version="seqyclean : $(seqyclean -h | grep Version)"

    seqyclean !{params.seqyclean_options} \
      -c !{params.seqyclean_contaminant_file} \
      -U !{reads} \
      -o seqyclean/!{sample}_cln \
      -gz \
      2>> $err_file >> $log_file
  '''
  } else if ( paired_single == 'paired' ) {
  '''
    mkdir -p seqyclean logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "seqyclean version: $(seqyclean -h | grep Version)" >> $log_file
    cleaner_version="seqyclean : $(seqyclean -h | grep Version)"

    seqyclean !{params.seqyclean_options} \
      -c !{params.seqyclean_contaminant_file} \
      -1 !{reads[0]} -2 !{reads[1]} \
      -o seqyclean/!{sample}_clean \
      -gz \
      2>> $err_file >> $log_file
  '''
  }
}
