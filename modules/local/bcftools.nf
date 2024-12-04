process BCFTOOLS {
  tag           "${meta.id}"
  publishDir    params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/bcftools:1.21'
  label         'process_single'

  when:
  params.bcftools_variants && (task.ext.when == null || task.ext.when)

  input:
  tuple val(meta), file(bam), file(reference_genome)

  output:
  tuple val(meta), file("bcftools_variants/*.vcf"), emit: vcf, optional: true
  path "bcftools_variants/*.vcf", emit: bcftools_variants_file, optional: true
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  shell:
  def args   = task.ext.args   ?: "-mv -Ov"
  def prefix = task.ext.prefix ?: "${meta.id}"
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
