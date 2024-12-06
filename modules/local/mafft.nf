process MAFFT {
  tag           "Multiple Sequence Alignment"
  label         "process_high"
  container     'staphb/mafft:7.526'

  input:
  file(fasta)
  file(reference_genome)

  output:
  path "mafft/mafft_aligned.fasta", emit: msa, optional: true
  path "mafft/*", emit: files
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: "${params.mafft_options}"
  def files  = fasta.join(" ")
  def prefix = task.ext.prefix ?: "mafft_aligned"
  def frags  = reference_genome ? "--addfragments" : ""
  def ref    = reference_genome ?: "" 
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
      ${frags} mafft/ultimate.fasta \
      ${ref} \
      > mafft/${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      mafft: \$(mafft --version 2>&1 | awk '{print \$1 "_" \$2}' | head -n 1)
      container: ${task.container}
    END_VERSIONS
  """
}
