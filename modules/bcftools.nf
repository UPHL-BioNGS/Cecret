process bcftools_variants {
  tag "${sample}"
  errorStrategy 'ignore'

  when:
  params.bcftools_variants

  input:
  tuple val(sample), file(bam), file(reference_genome)

  output:
  tuple val(sample), file("bcftools_variants/${sample}.vcf"),             emit: bcftools_variants_file
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p bcftools_variants logs/!{task.process}
    log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log
    bcftools --version >> $log

    bcftools mpileup -A -d !{params.mpileup_depth} -B -Q 0 -f !{reference_genome} !{bam} | \
      bcftools call -mv -Ov -o bcftools_variants/!{sample}.vcf | tee -a $log
  '''
}
