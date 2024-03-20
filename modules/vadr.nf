process vadr {
  tag           "QC metrics"
  label         "process_medium"
  publishDir    path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/vadr:1.6.3'

  //#UPHLICA maxForks 10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.vadr && (task.ext.when == null || task.ext.when)

  input:
  file(fasta)

  output:
  path "vadr/*", emit: vadr_files, optional: true
  path "vadr/vadr.vadr.sqa", emit: vadr_file,  optional: true
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log"
  path "versions.yml", emit: versions

  shell:
  def args      = task.ext.args      ?: "${params.vadr_options}"
  def trim_args = task.ext.trim_args ?: "${params.vadr_trim_options}"
  def fastas    = fasta.join(" ")
  """
    mkdir -p logs/${task.process}
    log=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    date > \$log
    v-annotate.pl -h | tee -a \$log

    for fasta in ${fastas}
    do
      lines=\$(grep -v ">" \$fasta | fold -w 75 | grep -e A -e G -e C -e T | wc -l | awk '{print \$1}' )
      if [ "\$lines" -gt 2 ] ; then cat \$fasta >> ultimate_fasta.fasta ; fi
    done

    if [ -f "ultimate_fasta.fasta" ]
    then
      fasta-trim-terminal-ambigs.pl ${trim_args} \
        ultimate_fasta.fasta > trimmed_ultimate.fasta
    fi

    if [ -s "trimmed_ultimate.fasta" ] &&  [ -f "trimmed_ultimate.fasta" ]
    then
      v-annotate.pl ${args} \
        --cpu ${task.cpus} \
        --noseqnamemax \
        --mkey ${params.vadr_reference} \
        --mdir ${params.vadr_mdir} \
        trimmed_ultimate.fasta \
        vadr \
        | tee -a \$log
    fi

    if [ -f "ultimate_fasta.fasta" ]   ; then cp ultimate_fasta.fasta vadr/combined.fasta  ; fi
    if [ -f "trimmed_ultimate.fasta" ] ; then cp trimmed_ultimate.fasta vadr/trimmed.fasta ; fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vadr: \$(v-annotate.pl -h | grep VADR | head -n 1 | awk '{print \$3 "_" \$4 "_" \$5}')
        container: ${task.container}
    END_VERSIONS
  """
}
