process bcftools_variants {
  tag "${sample}"

  when:
  params.bcftools_variants

  input:
  tuple val(sample), file(bam), file(reference_genome)

  output:
  tuple val(sample), file("bcftools_variants/${sample}.vcf"),             emit: bcftools_variants_file
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}",  emit: log
  tuple val(sample), env(variants_num),                                   emit: bcftools_variants_results

  shell:
  '''
    mkdir -p bcftools_variants logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    bcftools --version >> $log_file

    bcftools mpileup -A -d !{params.mpileup_depth} -B -Q 0 -f !{reference_genome} !{bam} 2>> $err_file | \
      bcftools call -mv -Ov -o bcftools_variants/!{sample}.vcf 2>> $err_file >> $log_file

    variants_num=$(grep -v "#" bcftools_variants/!{sample}.vcf | wc -l)
    if [ -z "$variants_num" ] ; then variants_num="0" ; fi
  '''
}
