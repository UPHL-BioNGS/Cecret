process fastqc {
  tag        "${sample}"
  label      "process_single"
  publishDir "${params.outdir}", mode: 'copy'
  container  'staphb/fastqc:0.12.1'

  //#UPHLICA maxForks 10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '45m'

  when:
  params.fastqc && sample != null

  input:
  tuple val(sample), file(fastq), val(type)

  output:
  path "fastqc/*.html",                   emit: files
  path "fastqc/*_fastqc.zip",             emit: fastqc_files
  path "fastqc/${sample}_fastq_name.csv", emit: fastq_name
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p fastqc logs/!{task.process}
    log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log
    fastqc --version >> $log

    fastqc !{params.fastqc_options} \
      --outdir fastqc \
      --threads !{task.cpus} \
      !{fastq} \
      | tee -a $log

    echo "sample,fastq_1,fastq_2"             > fastqc/!{sample}_fastq_name.csv
    echo "!{sample},!{fastq[0]},!{fastq[1]}" >> fastqc/!{sample}_fastq_name.csv
  '''
}
