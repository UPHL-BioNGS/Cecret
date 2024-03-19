process freyja_variants {
  tag           "${sample}"
  label         "process_medium"
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  publishDir    path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/freyja:1.4.9-03_18_2024-00-45-2024-03-19'

  //#UPHLICA maxForks 10
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.freyja && (task.ext.when == null || task.ext.when)

  input:
  tuple val(sample), file(bam), file(reference_genome)

  output:
  tuple val(sample), file("freyja/${sample}_{depths,variants}.tsv"), optional: true, emit: variants
  path "freyja/${sample}*", optional: true, emit: files
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"
  path "versions.yml", emit: versions

  shell:
  def args   = task.ext.args   ?: "${params.freyja_variants_options}"
  def prefix = task.ext.prefix ?: "${sample}"
  """
    mkdir -p freyja logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    date > \$log
    freyja --version >> \$log

    freyja variants ${args} \
      ${bam} \
      --variants freyja/${prefix}_variants.tsv \
      --depths freyja/${prefix}_depths.tsv \
      --ref ${reference_genome} \
      | tee -a \$log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      freyja: \$(freyja --version | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}

process freyja_demix {
  tag           "${sample}"
  label         "process_medium"
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  publishDir    path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/freyja:1.4.9-03_18_2024-00-45-2024-03-19'


  //#UPHLICA maxForks 10
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.freyja && (task.ext.when == null || task.ext.when)

  input:
  tuple val(sample), file(variants)

  output:
  path "freyja/${sample}_demix.tsv", optional: true, emit: demix
  path "freyja/${sample}*",          optional: true, emit: files
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"
  path "versions.yml", emit: versions

  shell:
  def args   = task.ext.args   ?: "${params.freyja_demix_options}"
  def prefix = task.ext.prefix ?: "${sample}"
  """
    mkdir -p freyja logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    date > \$log
    freyja --version >> \$log

    freyja demix ${args} \
      ${variants[1]} \
      ${variants[0]} \
      --output freyja/${prefix}_demix.tsv \
      | tee -a \$log

    freyja --help

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      freyja: \$(freyja --version | awk '{print \$NF}')
      barcode: \$(freyja demix --version | head -n 2 | tail -n 1 )
      container: ${task.container}
    END_VERSIONS
  """
}

process freyja_aggregate {
  tag        "Aggregating results from freyja"
  label      "process_single"
  publishDir path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container  'staphb/freyja:1.4.9-03_18_2024-00-45-2024-03-19'

  //#UPHLICA maxForks 10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.freyja_aggregate && (task.ext.when == null || task.ext.when)

  input:
  file(demix)
  file(graphs)

  output:
  path "freyja/*", emit: files
  path "freyja/aggregated-freyja.tsv", emit: aggregated_freyja_file, optional: true
  path "freyja/*mqc.png", emit: for_multiqc
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log"
  path "versions.yml", emit: versions

  shell:
  def args   = task.ext.args   ?: "${params.freyja_aggregate_options}"
  def pltarg = task.ext.pltarg ?: "${params.freyja_plot_options}"
  def filtyp = task.ext.filtyp ?: "${params.freyja_plot_filetype}"
  def prefix = task.ext.prefix ?: "aggregated-freyja"
  def files  = demix.join(" ")
  """
    mkdir -p freyja logs/${task.process} tmp
    log=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    date > \$log

    mv ${files} tmp/.

    freyja aggregate ${args} \
      tmp/ \
      --output freyja/${prefix}.tsv \
      | tee -a \$log

    freyja plot ${pltarg} \
      freyja/${prefix}.tsv \
      --output freyja/${prefix}.${filtyp} \
      | tee -a \$log

    python3 ${graphs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      freyja: \$(freyja --version | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}
