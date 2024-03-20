process seqyclean {
  tag           "${sample}"
  label         "process_single"
  publishDir    path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/seqyclean:1.10.09'

  //#UPHLICA maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '45m'

  when:
  sample != null && (task.ext.when == null || task.ext.when)

  input:
  tuple val(sample), file(reads), val(paired_single)

  output:
  tuple val(sample), file("seqyclean/${sample}_{cln_SE,clean_PE1,clean_PE2}.fastq.gz"), optional: true, emit: clean_reads
  path "seqyclean/${sample}_clean_SummaryStatistics.tsv", optional: true, emit: seqyclean_files_collect_paired
  path "seqyclean/${sample}_cln_SummaryStatistics.tsv", optional: true, emit: seqyclean_files_collect_single
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"
  tuple val("${params.cleaner}"), env(cleaner_version), emit: cleaner_version
  path "versions.yml", emit: versions

  shell:
  def args   = task.ext.args   ?: "${params.seqyclean_options}"
  def prefix = task.ext.prefix ?: "${sample}"
  if ( paired_single == "single" ) {
  """
    mkdir -p seqyclean logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log
    echo "seqyclean version: \$(seqyclean -h | grep Version)" >> \$log
    cleaner_version="seqyclean : \$(seqyclean -h | grep Version)"

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
  } else if ( paired_single == 'paired' ) {
  """
    mkdir -p seqyclean logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log
    echo "seqyclean version: \$(seqyclean -h | grep Version)" >> \$log
    cleaner_version="seqyclean : \$(seqyclean -h | grep Version)"

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
