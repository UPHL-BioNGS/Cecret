process SEQYCLEAN {
  tag           "${meta.id}"
  label         "process_single"
  publishDir    path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/seqyclean:1.10.09'

  when:
  task.ext.when == null || task.ext.when

  input:
  tuple val(meta), file(reads)
  
  output:
  tuple val(meta), file("seqyclean/*_{cln_SE,clean_PE1,clean_PE2}.fastq.gz"), optional: true, emit: clean_reads
  path "seqyclean/*_clean_SummaryStatistics.tsv", optional: true, emit: seqyclean_files_collect_paired
  path "seqyclean/*_cln_SummaryStatistics.tsv", optional: true, emit: seqyclean_files_collect_single
  path "logs/${task.process}/*.log", emit: log
  tuple val("${params.cleaner}"), env(cleaner_version), emit: cleaner_version
  path "versions.yml", emit: versions

  shell:
  def args   = task.ext.args   ?: "${params.seqyclean_options}"
  def prefix = task.ext.prefix ?: "${meta.id}"
  if ( meta.single_end ) {
  """
    mkdir -p seqyclean logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log
    echo "seqyclean version: \$(seqyclean -h | grep Version)" >> \$log
    cleaner_version=\$(seqyclean -h | grep Version | head -n 1 | awk '{print \$2 "_" \$3}' )

    seqyclean ${args} \
      -c ${params.seqyclean_contaminant_file} \
      -U ${reads} \
      -o seqyclean/${prefix}_cln \
      -gz \
      | tee -a \$log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      seqyclean: \$(seqyclean -h | grep Version | head -n 1 | awk '{print \$2 "_" \$3}' )
      container: ${task.container}
    END_VERSIONS
  """
  } else {
  """
    mkdir -p seqyclean logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log
    echo "seqyclean version: \$(seqyclean -h | grep Version)" >> \$log
    cleaner_version=\$(seqyclean -h | grep Version | head -n 1 | awk '{print \$2 "_" \$3}' )

    seqyclean ${args} \
      -c ${params.seqyclean_contaminant_file} \
      -1 ${reads[0]} -2 ${reads[1]} \
      -o seqyclean/${prefix}_clean \
      -gz \
      | tee -a \$log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      seqyclean: \$(seqyclean -h | grep Version | head -n 1 | awk '{print \$2 "_" \$3}' )
      container: ${task.container}
    END_VERSIONS
  """
  }
}
