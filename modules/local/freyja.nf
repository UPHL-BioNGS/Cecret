process FREYJA {
  tag           "${meta.id}"
  label         "process_medium"
  container     'staphb/freyja:1.5.3-06_23_2025-00-42-2025-06-23'

  input:
  tuple val(meta), file(bam), file(reference_genome)

  output:
  tuple val(meta), file("freyja/*_{depths,variants}.tsv"), optional: true, emit: variants
  path("freyja/*_demix.tsv"), optional: true, emit: demix
  path "freyja/*", optional: true, emit: files
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: "${params.freyja_variants_options}"
  def prefix = task.ext.prefix ?: "${meta.id}"
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

    freyja demix ${args} \
      freyja/${prefix}_variants.tsv \
      freyja/${prefix}_depths.tsv \
      --output freyja/${prefix}_demix.tsv \
      | tee -a \$log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      freyja: \$(freyja --version | awk '{print \$NF}')
      barcode: \$(freyja demix --version | head -n 2 | tail -n 1 )
      container: ${task.container}
    END_VERSIONS
  """
}

process FREYJA_AGGREGATE {
  tag        "Aggregating results from freyja"
  label      "process_single"
  container  'staphb/freyja:1.5.3-06_23_2025-00-42-2025-06-23'

  input:
  file(demix)
  file(graphs)

  output:
  path "freyja/*", emit: files
  path "freyja/aggregated-freyja.tsv", emit: aggregated_freyja_file, optional: true
  path "freyja/*mqc.png", emit: for_multiqc
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
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

    python3 freyja_graphs.py

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      freyja: \$(freyja --version | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}
