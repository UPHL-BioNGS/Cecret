process nextalign {
  tag        "Multiple Sequence Alignment"
  label      "process_high"
  publishDir "${params.outdir}", mode: 'copy'
  container  'nextstrain/nextalign:2.14.0'

  //#UPHLICA maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  input:
  file(consensus)
  path(dataset)

  output:
  path "nextalign/nextalign.aligned.fasta",                                     emit: msa
  path "nextalign/{*.fasta,nextalign.*.csv}",                                   emit: files
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p nextalign logs/!{task.process}
    log=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log

    date > $log
    echo "nextalign version:" >> $log
    nextalign --version >> $log

    for fasta in !{consensus}
    do
      cat $fasta >> nextalign/ultimate.fasta
    done

    nextalign run !{params.nextalign_options} \
      --input-ref=!{dataset}/reference.fasta \
      --genemap=!{dataset}/genemap.gff \
      --jobs !{task.cpus} \
      --output-all=nextalign/ \
      nextalign/ultimate.fasta \
      | tee -a $log
  '''
}
