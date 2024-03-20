process bcftools_variants {
  tag           "${sample}"
  publishDir    params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  container     'staphb/bcftools:1.19'
  label         'process_single'

  //#UPHLICA maxForks 10
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '45m'

  when:
  params.bcftools_variants && (task.ext.when == null || task.ext.when)

  input:
  tuple val(sample), file(bam), file(reference_genome)

  output:
  tuple val(sample), file("bcftools_variants/${sample}.vcf"), emit: vcf, optional: true
  path "bcftools_variants/${sample}.vcf", emit: bcftools_variants_file, optional: true
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"
  path "versions.yml", emit: versions

  shell:
  def args   = task.ext.args   ?: "-mv -Ov"
  def prefix = task.ext.prefix ?: "${sample}"
  """
    mkdir -p bcftools_variants logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log
    bcftools --version >> \$log

    bcftools mpileup -A -d ${params.mpileup_depth} -B -Q 0 -f ${reference_genome} ${bam} | \
      bcftools call ${args} -o bcftools_variants/${prefix}.vcf | tee -a \$log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      bcftools: \$(bcftools --version | head -n 1 | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}
