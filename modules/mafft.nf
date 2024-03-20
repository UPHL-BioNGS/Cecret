process mafft {
  tag           "Multiple Sequence Alignment"
  label         "process_high"
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  publishDir    path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/mafft:7.505'

  //#UPHLICA maxForks 10
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '2h'

  when:
  task.ext.when == null || task.ext.when

  input:
  file(fasta)
  file(reference_genome)

  output:
  path "mafft/mafft_aligned.fasta", emit: msa, optional: true
  path "mafft/*", emit: files
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log"
  path "versions.yml", emit: versions

  shell:
  def args   = task.ext.args   ?: "${params.mafft_options}"
  def files  = fasta.join(" ")
  def prefix = task.ext.prefix ?: "mafft_aligned"
  """
    mkdir -p mafft logs/${task.process}
    log=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    date > \$log
    echo "mafft version:" >> \$log
    mafft --version 2>&1 >> \$log

    for fasta in ${files}
    do
      cat \$fasta >> mafft/ultimate.fasta
    done

    mafft --auto \
      ${args} \
      --thread ${task.cpus} \
      --addfragments mafft/ultimate.fasta \
      ${reference_genome} \
      > mafft/${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      mafft: \$(mafft --version 2>&1 | awk '{print \$1 "_" \$2}' | head -n 1)
      container: ${task.container}
    END_VERSIONS
  """
}
