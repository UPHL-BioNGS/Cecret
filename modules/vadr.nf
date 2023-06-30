process vadr {
  tag        "QC metrics"
  label      "medcpus"
  publishDir "${params.outdir}", mode: 'copy'
  container  'staphb/vadr:1.5.1'

  //#UPHLICA maxForks 10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.vadr

  input:
  file(fasta)

  output:
  path "vadr/*",              emit: vadr_files
  path "vadr/vadr.vadr.sqa",  emit: vadr_file
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p logs/!{task.process}
    log=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log

    date > $log
    vadr --version | tee -a $log
    v-annotate.pl -h | tee -a $log

    for fasta in !{fasta}
    do
      lines=$(grep -v ">" $fasta | fold -w 75 | grep -e A -e G -e C -e T | wc -l | awk '{print $1}' )
      if [ "$lines" -gt 2 ] ; then cat $fasta >> ultimate_fasta.fasta ; fi
    done

    fasta-trim-terminal-ambigs.pl !{params.vadr_trim_options} \
      ultimate_fasta.fasta > trimmed_ultimate_fasta.fasta

    if [ -s "trimmed_ultimate_fasta.fasta" ]
    then
      v-annotate.pl !{params.vadr_options} \
        --cpu !{task.cpus} \
        --noseqnamemax \
        --mkey !{params.vadr_reference} \
        --mdir !{params.vadr_mdir} \
        trimmed_ultimate_fasta.fasta \
        vadr \
        | tee -a $log
    fi
    cp ultimate_fasta.fasta vadr/combined.fasta
    cp trimmed_ultimate_fasta.fasta vadr/trimmed.fasta
  '''
}
