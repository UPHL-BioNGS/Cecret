process bcftools_variants {
  tag           "${sample}"
  publishDir    "${params.outdir}", mode: 'copy'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  container     'staphb/bcftools:1.18'
  label         'process_single'

  //#UPHLICA maxForks 10
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '45m'

  when:
  params.bcftools_variants

  input:
  tuple val(sample), file(bam), file(reference_genome)

  output:
  tuple val(sample), file(bam), file(reference_genome), file("bcftools_variants/${sample}.vcf"), emit: vcf
  path "bcftools_variants/${sample}.vcf", emit: bcftools_variants_file
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
