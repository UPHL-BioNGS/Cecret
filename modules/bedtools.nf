process bedtools_multicov {
  tag "${sample}"

  when:
  params.bedtools_multicov

  input:
  tuple val(sample), file(bam), file(bai), file(amplicon_bed)

  output:
  path "multicov/${sample}.multicov.txt",                                 emit: multicov
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p multicov logs/!{task.process}
    log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    date > $log
    bedtools --version >> $log

    bedtools multicov !{params.bedtools_multicov_options} \
      -bams !{bam} \
      -bed !{amplicon_bed} \
      >> multicov/!{sample}.multicov.txt
  '''
}
