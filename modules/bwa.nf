process bwa {
  tag         "${sample}"
  label       "process_high"
  publishDir  path: "${params.outdir}", mode: 'copy', pattern: 'logs/*/*log'
  container   'staphb/bwa:0.7.17'
  
  //#UPHLICA maxForks 10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  input:
  tuple val(sample), file(reads), file(reference_genome)

  output:
  tuple val(sample), file("bwa/${sample}.sam"),                     emit: sam
  tuple val("${params.aligner}"), env(bwa_version),                 emit: aligner_version
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p bwa logs/!{task.process}
    log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log
    echo "bwa $(bwa 2>&1 | grep Version )" >> $log
    bwa_version="bwa : "$(bwa 2>&1 | grep Version)

    # index the reference fasta file
    bwa index !{reference_genome}

    # bwa mem command
    bwa mem -t !{task.cpus} !{reference_genome} !{reads} > bwa/!{sample}.sam
  '''
}
