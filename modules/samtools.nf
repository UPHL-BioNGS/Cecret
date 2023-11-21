process samtools_stats {
  tag        "${sample}"
  label      "process_single"
  publishDir "${params.outdir}", mode: 'copy'
  container  'staphb/samtools:1.18'
  
  //#UPHLICA maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.samtools_stats

  input:
  tuple val(sample), file(bam)

  output:
  path "samtools_stats/${sample}.stats.txt", emit: samtools_stats_files
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p samtools_stats logs/!{task.process}
    log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    date > $log
    samtools --version >> $log

    samtools stats !{params.samtools_stats_options} !{bam} > samtools_stats/!{sample}.stats.txt
  '''
}

process samtools_coverage {
  tag        "${sample}"
  label      "process_single"
  publishDir "${params.outdir}", mode: 'copy'
  container  'staphb/samtools:1.18'
  
  //#UPHLICA maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.samtools_coverage

  input:
  tuple val(sample), file(bam)

  output:
  path "samtools_coverage/${sample}.cov.{txt,hist}", emit: files
  path "samtools_coverage/${sample}.cov.txt", emit: samtools_coverage
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p samtools_coverage logs/!{task.process}
    log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    date > $log
    samtools --version >> $log

    samtools coverage !{params.samtools_coverage_options} !{bam} -m -o samtools_coverage/!{sample}.cov.hist | tee -a $log
    samtools coverage !{params.samtools_coverage_options} !{bam} | awk -v sample=!{sample} '{print sample "\\t" $0 }' | sed '0,/!{sample}/s//sample/' > samtools_coverage/!{sample}.cov.txt  | tee -a $log
  '''
}

process samtools_flagstat {
  tag        "${sample}"
  label      "process_single"
  publishDir "${params.outdir}", mode: 'copy'
  container  'staphb/samtools:1.18'
  
  //#UPHLICA maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  input:
  tuple val(sample), file(bam)

  when:
  params.samtools_flagstat

  output:
  path "samtools_flagstat/${sample}.flagstat.txt", emit: samtools_flagstat_files
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p samtools_flagstat logs/!{task.process}
    log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    date > $log
    samtools --version >> $log

    samtools flagstat !{params.samtools_flagstat_options} \
      !{bam} | \
      tee samtools_flagstat/!{sample}.flagstat.txt
  '''
}

process samtools_depth {
  tag        "${sample}"
  label      "process_single"
  publishDir "${params.outdir}", mode: 'copy'
  container  'staphb/samtools:1.18'
  
  //#UPHLICA maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  input:
  tuple val(sample), file(bam)

  when:
  params.samtools_depth

  output:
  path "samtools_depth/${sample}.depth.txt", emit: file
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p samtools_depth logs/!{task.process}
    log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    date > $log
    samtools --version >> $log

    samtools depth !{params.samtools_depth_options} \
      !{bam} > samtools_depth/!{sample}.depth.txt
  '''
}

process samtools_ampliconstats {
  tag        "${sample}"
  label      "process_single"
  publishDir "${params.outdir}", mode: 'copy'
  container  'staphb/samtools:1.18'
  
  //#UPHLICA maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.samtools_ampliconstats && ( params.trimmer != 'none' )

  input:
  tuple val(sample), file(bam), file(primer_bed)

  output:
  tuple val(sample), file("samtools_ampliconstats/${sample}_ampliconstats.txt"), emit: samtools_ampliconstats_files
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p samtools_ampliconstats logs/!{task.process}
    log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    date > $log
    samtools --version >> $log

    samtools ampliconstats !{params.samtools_ampliconstats_options} \
      !{primer_bed} \
      !{bam} > samtools_ampliconstats/!{sample}_ampliconstats.txt
  '''
}

process samtools_plot_ampliconstats {
  tag           "${sample}"
  label      "process_single"
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  publishDir    "${params.outdir}", mode: 'copy'
  container     'staphb/samtools:1.18'

  //#UPHLICA maxForks 10
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.samtools_plot_ampliconstats

  input:
  tuple val(sample), file(ampliconstats)

  output:
  path "samtools_plot_ampliconstats/${sample}*", emit: files
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p samtools_plot_ampliconstats/!{sample} logs/!{task.process}
    log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    date > $log
    samtools --version >> $log

    plot-ampliconstats !{params.samtools_plot_ampliconstats_options} \
      samtools_plot_ampliconstats/!{sample} \
      !{ampliconstats}
  '''
}

process samtools_sort {
  tag        "${sample}"
  label      "process_high"
  publishDir "${params.outdir}", mode: 'copy'
  container  'staphb/samtools:1.18'

  //#UPHLICA maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  input:
  tuple val(sample), file(sam)

  output:
  tuple val(sample), file("aligned/${sample}.sorted.bam"), file("aligned/${sample}.sorted.bam.bai"),  emit: bam_bai
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p aligned logs/!{task.process}
    log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log
    samtools --version >> $log

    samtools sort -@ !{task.cpus} !{sam} | \
      samtools view -F 4 -o aligned/!{sample}.sorted.bam | tee -a $log

    # indexing the bams
    samtools index aligned/!{sample}.sorted.bam | tee -a $log
  '''
}

process samtools_filter {
  tag        "${sample}"
  label      "process_single"
  publishDir "${params.outdir}", mode: 'copy'
  container  'staphb/samtools:1.18'

  //#UPHLICA maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.filter

  input:
  tuple val(sample), file(sam)

  output:
  tuple val(sample), file("filter/${sample}_filtered_{R1,R2}.fastq.gz"), optional: true, emit: filtered_reads
  path "filter/${sample}_filtered_unpaired.fastq.gz",                    optional: true, emit: unpaired
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p filter logs/!{task.process}
    log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log
    samtools --version >> $log

    samtools sort -n !{sam} | \
      samtools fastq -F 4 !{params.filter_options} \
      -s filter/!{sample}_filtered_unpaired.fastq.gz \
      -1 filter/!{sample}_filtered_R1.fastq.gz \
      -2 filter/!{sample}_filtered_R2.fastq.gz \
      | tee -a $log
  '''
}

process samtools_ampliconclip {
  tag        "${sample}"
  label      "process_single"
  publishDir "${params.outdir}", mode: 'copy'
  container  'staphb/samtools:1.18'
  
  //#UPHLICA maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  input:
  tuple val(sample), file(bam), file(primer_bed)

  output:
  tuple val(sample), file("ampliconclip/${sample}.primertrim.sorted.bam"), file("ampliconclip/${sample}.primertrim.sorted.bam.bai"), emit: bam_bai
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"                                                         
  tuple val("samtools ampliconclip"), env(trimmer_version),                                                                          emit: trimmer_version

  shell:
  '''
    mkdir -p ampliconclip logs/!{task.process}
    log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log
    samtools --version >> $log
    trimmer_version="samtools ampliconclip : $(samtools --version | head -n 1)"

    # trimming the reads
    samtools ampliconclip !{params.samtools_ampliconclip_options} -b !{primer_bed} !{bam} | \
      samtools sort |  \
      samtools view -F 4 -o ampliconclip/!{sample}.primertrim.sorted.bam | tee -a $log

    samtools index ampliconclip/!{sample}.primertrim.sorted.bam | tee -a $log
  '''
}

process samtools_markdup {
  tag        "${sample}"
  label      "process_single"
  publishDir "${params.outdir}", mode: 'copy'
  container  'staphb/samtools:1.18'
  
  //#UPHLICA maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  input:
  tuple val(sample), val(type), file(sam) 

  output:
  tuple val(sample), file("markdup/${sample}.markdup.sorted.bam"), file("markdup/${sample}.markdup.sorted.bam.bai"),  emit: bam_bai
  path "markdup/${sample}_markdupstats.txt",                                                                          emit: stats
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"
  
  shell:
  if ( type == 'single' ) {
  '''
    mkdir -p markdup logs/!{task.process}
    log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log
    samtools --version >> $log

    samtools sort !{sam} | \
      samtools markdup !{params.samtools_markdup_options} -@ !{task.cpus} -s -f markdup/!{sample}_markdupstats.txt - markdup/!{sample}.markdup.sorted.bam | tee -a $log

    samtools index markdup/!{sample}.markdup.sorted.bam | tee -a $log
  '''
  } else if (type == 'paired') {
  '''
    mkdir -p markdup logs/!{task.process}
    log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log
    samtools --version >> $log

    samtools sort -n !{sam} | \
      samtools fixmate !{params.samtools_fixmate_options} -m -@ !{task.cpus} - - | \
      samtools sort | \
      samtools markdup !{params.samtools_markdup_options} -@ !{task.cpus} -s -f markdup/!{sample}_markdupstats.txt - markdup/!{sample}.markdup.sorted.bam | tee -a $log

    samtools index markdup/!{sample}.markdup.sorted.bam | tee -a $log
  '''
  }
}

