process samtools_stats {
  tag "${sample}"

  when:
  params.samtools_stats

  input:
  tuple val(sample), file(bam)

  output:
  path "samtools_stats/${sample}.stats.txt",                              emit: samtools_stats_files
  tuple val(sample), env(insert_size_after_trimming),                     emit: samtools_stats_after_size_results
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}",  emit: log

  shell:
  '''
    mkdir -p samtools_stats logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    samtools stats !{params.samtools_stats_options} !{bam} 2>> $err_file > samtools_stats/!{sample}.stats.txt

    insert_size_after_trimming=$(grep "insert size average" samtools_stats/!{sample}.stats.txt | cut -f 3)
    if [ -z "$insert_size_after_trimming" ] ; then insert_size_after_trimming=0 ; fi
  '''
}

process samtools_coverage {
  tag "${sample}"

  when:
  params.samtools_coverage

  input:
  tuple val(sample), file(bam)

  output:
  path "samtools_coverage/${sample}.cov.{txt,hist}",                      emit: files
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}",  emit: log
  tuple val(sample), env(coverage),                                       emit: samtools_coverage_results
  tuple val(sample), env(depth),                                          emit: samtools_covdepth_results

  shell:
  '''
    mkdir -p samtools_coverage logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    samtools coverage !{params.samtools_coverage_options} !{bam} -m -o samtools_coverage/!{sample}.cov.hist 2>> $err_file >> $log_file
    samtools coverage !{params.samtools_coverage_options} !{bam}    -o samtools_coverage/!{sample}.cov.txt  2>> $err_file >> $log_file

    coverage=$(cut -f 6 samtools_coverage/!{sample}.cov.txt | tail -n 1)
    depth=$(cut    -f 7 samtools_coverage/!{sample}.cov.txt | tail -n 1)

    if [ -z "$coverage" ] ; then coverage="0" ; fi
    if [ -z "$depth" ] ; then depth="0" ; fi
  '''
}

process samtools_flagstat {
  tag "${sample}"

  input:
  tuple val(sample), file(bam)

  when:
  params.samtools_flagstat

  output:
  path "samtools_flagstat/${sample}.flagstat.txt",                        emit: samtools_flagstat_files
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}",  emit: log

  shell:
  '''
    mkdir -p samtools_flagstat logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    samtools flagstat !{params.samtools_flagstat_options} \
      !{bam} \
      2>> $err_file | tee samtools_flagstat/!{sample}.flagstat.txt
  '''
}

process samtools_depth {
  tag "${sample}"

  input:
  tuple val(sample), file(bam)

  when:
  params.samtools_depth

  output:
  path "samtools_depth/${sample}.depth.txt",                              emit: file
  tuple val(sample), env(depth),                                          emit: samtools_depth_results
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}",  emit: log

  shell:
  '''
    mkdir -p samtools_depth logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    samtools depth !{params.samtools_depth_options} \
      !{bam} \
      2>> $err_file > samtools_depth/!{sample}.depth.txt

    depth=$(awk '{ if ($3 > !{params.minimum_depth} ) print $0 }' samtools_depth/!{sample}.depth.txt | wc -l )
    if [ -z "$depth" ] ; then depth="0" ; fi
  '''
}

process samtools_ampliconstats {
  tag "${sample}"

 when:
 params.samtools_ampliconstats && ( params.trimmer != 'none' )

  input:
  tuple val(sample), file(bam), file(primer_bed)

  output:
  tuple val(sample), file("samtools_ampliconstats/${sample}_ampliconstats.txt"), emit: samtools_ampliconstats_files
  tuple val(sample), env(num_failed_amplicons),                                  emit: samtools_ampliconstats_results
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}",         emit: log

  shell:
  '''
    mkdir -p samtools_ampliconstats logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    samtools ampliconstats !{params.samtools_ampliconstats_options} \
      !{primer_bed} \
      !{bam} \
      2>> $err_file > samtools_ampliconstats/!{sample}_ampliconstats.txt

    num_failed_amplicons=$(grep ^FREADS samtools_ampliconstats/!{sample}_ampliconstats.txt | cut -f 2- | tr '\t' '\n' | awk '{ if ($1 < 20) print $0 }' | wc -l)
    if [ -z "$num_failed_amplicons" ] ; then num_failed_amplicons=0 ; fi
  '''
}

process samtools_plot_ampliconstats {
  tag "${sample}"
  errorStrategy 'ignore'

  when:
  params.samtools_plot_ampliconstats

  input:
  tuple val(sample), file(ampliconstats)

  output:
  path "samtools_plot_ampliconstats/${sample}*",                          emit: files
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}" , emit: log

  shell:
  '''
    mkdir -p samtools_plot_ampliconstats/!{sample} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    plot-ampliconstats !{params.samtools_plot_ampliconstats_options} \
      samtools_plot_ampliconstats/!{sample} \
      !{ampliconstats}
  '''
}

process samtools_sort {
  tag "${sample}"
  label "maxcpus"

  input:
  tuple val(sample), file(sam)

  output:
  tuple val(sample), file("aligned/${sample}.sorted.bam"),                                            emit: bam
  tuple val(sample), file("aligned/${sample}.sorted.bam"), file("aligned/${sample}.sorted.bam.bai"),  emit: bam_bai
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}",                              emit: log

  shell:
  '''
    mkdir -p aligned logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    # NOTE: the -@ flag is used for thread count
    samtools sort -@ !{task.cpus} !{sam} 2>> $err_file | \
      samtools view -F 4 -o aligned/!{sample}.sorted.bam 2>> $err_file >> $log_file

    # indexing the bams
    samtools index aligned/!{sample}.sorted.bam 2>> $err_file >> $log_file
  '''
}

process samtools_filter {
  tag "${sample}"

  when:
  params.filter

  input:
  tuple val(sample), file(sam)

  output:
  tuple val(sample), file("filter/${sample}_filtered_{R1,R2}.fastq.gz"), optional: true, emit: filtered_reads
  path "filter/${sample}_filtered_unpaired.fastq.gz",                    optional: true, emit: unpaired
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}",                 emit: log

  shell:
  '''
    mkdir -p filter logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    samtools sort -n !{sam} 2>> $err_file | \
      samtools fastq -F 4 !{params.filter_options} \
      -s filter/!{sample}_filtered_unpaired.fastq.gz \
      -1 filter/!{sample}_filtered_R1.fastq.gz \
      -2 filter/!{sample}_filtered_R2.fastq.gz \
      2>> $err_file >> $log_file
  '''
}

process samtools_ampliconclip {
  tag "${sample}"

  input:
  tuple val(sample), file(bam), file(primer_bed)

  output:
  tuple val(sample), file("ampliconclip/${sample}.primertrim.sorted.bam"),                                                           emit: trimmed_bam
  tuple val(sample), file("ampliconclip/${sample}.primertrim.sorted.bam"), file("ampliconclip/${sample}.primertrim.sorted.bam.bai"), emit: bam_bai
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}",                                                             emit: log
  tuple val(sample), env(trimmer_version),                                                                                           emit: trimmer_version

  shell:
  '''
    mkdir -p ampliconclip logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file
    trimmer_version="samtools ampliconclip : $(samtools --version | head -n 1)"

    # trimming the reads
    samtools ampliconclip !{params.samtools_ampliconclip_options} -b !{primer_bed} !{bam} 2>> $err_file | \
      samtools sort 2>> $err_file |  \
      samtools view -F 4 -o ampliconclip/!{sample}.primertrim.sorted.bam 2>> $err_file >> $log_file

    samtools index ampliconclip/!{sample}.primertrim.sorted.bam 2>> $err_file >> $log_file
  '''
}

process samtools_markdup {
  tag "${sample}"

  input:
  tuple val(sample), val(type), file(sam) 

  output:
  tuple val(sample), file("markdup/${sample}.markdup.sorted.bam"),                                                    emit: bam
  tuple val(sample), file("markdup/${sample}.markdup.sorted.bam"), file("markdup/${sample}.markdup.sorted.bam.bai"),  emit: bam_bai
  path "markdup/${sample}_markdupstats.txt",                                                                          emit: stats
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}",                                              emit: log
  
  shell:
  if ( type == 'single' ) {
  '''
    mkdir -p markdup logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    samtools sort !{sam} 2>> $err_file | \
      samtools markdup !{params.samtools_markdup_options} -@ !{task.cpus} -s -f markdup/!{sample}_markdupstats.txt - markdup/!{sample}.markdup.sorted.bam 2>> $err_file >> $log_file

    samtools index markdup/!{sample}.markdup.sorted.bam 2>> $err_file >> $log_file
  '''
  } else if (type == 'paired') {
  '''
    mkdir -p markdup logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    samtools sort -n !{sam} 2>> $err_file | \
      samtools fixmate !{params.samtools_fixmate_options} -m -@ !{task.cpus} - - 2>> $err_file | \
      samtools sort 2>> $err_file | \
      samtools markdup !{params.samtools_markdup_options} -@ !{task.cpus} -s -f markdup/!{sample}_markdupstats.txt - markdup/!{sample}.markdup.sorted.bam 2>> $err_file >> $log_file

    samtools index markdup/!{sample}.markdup.sorted.bam 2>> $err_file >> $log_file
  '''
  }
}