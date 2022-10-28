process bedtools_multicov {
  tag "${sample}"

  when:
  params.bedtools_multicov

  input:
  tuple val(sample), file(bam), file(bai), file(amplicon_bed)

  output:
  path "multicov/${sample}.multicov.txt",                                 emit: multicov
  tuple val(sample), env(num_failed_amplicons),                           emit: bedtools_results
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

    result_column=$(head -n 1 multicov/!{sample}.multicov.txt | awk '{print NF}' )
    num_failed_amplicons=$(cat multicov/!{sample}.multicov.txt | tr ' ' '\t' | cut -f $result_column | awk '{ if ( $1 < 20 ) print $0 }' | wc -l )
    if [ -z "$num_failed_amplicons" ] ; then num_failed_amplicons="NA" ; fi
  '''
}
