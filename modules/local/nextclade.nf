process NEXTCLADE_DATASET {
  tag        "Downloading Dataset"
  label      "process_medium"
  container  'nextstrain/nextclade:3.10.1'

  output:
  path "dataset", emit: dataset
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args ?: "${params.nextclade_dataset}"
  """
    mkdir -p nextclade dataset logs/${task.process}
    log=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    date > \$log
    nextclade --version >> \$log
    nextclade_version=\$(nextclade --version)

    echo "Getting nextclade dataset for ${args}" | tee -a \$log
    nextclade dataset list | tee -a \$log

    nextclade dataset get --name ${args} --output-dir dataset

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      nextclade: \$(nextclade --version | awk '{print \$NF}')
      tag: \$(grep "tag" dataset/pathogen.json | grep tag | sed 's/\"//g' | sed 's/,//g' | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}

process NEXTCLADE {
  tag        "Clade Determination"
  label      "process_medium"
  container  'nextstrain/nextclade:3.10.1'

  input:
  file(fasta)
  path(dataset)

  output:
  path "nextclade/nextclade.csv", emit: nextclade_file
  path "nextclade/*", emit: results
  tuple file("nextclade/nextclade.aligned.fasta"), file("nextclade/nextclade.nwk"), emit: prealigned, optional: true
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: "${params.nextclade_options}"
  def files  = fasta.join(" ")
  def prefix = task.ext.prefix ?: "combined"
  """
    mkdir -p nextclade dataset logs/${task.process}
    log=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    date > \$log
    nextclade --version >> \$log
    nextclade_version=\$(nextclade --version)

    for fasta in ${files}
    do
      cat \$fasta >> ultimate_fasta.fasta
    done

    nextclade run ${args} \
      --input-dataset ${dataset} \
      --output-all=nextclade/ \
      --jobs ${task.cpus} \
      ultimate_fasta.fasta \
      | tee -a \$log

    cp ultimate_fasta.fasta nextclade/${prefix}.fasta

    if [ -f "dataset/pathogen.json" ]
    then
      tag=\$(grep "tag" dataset/pathogen.json | grep tag | sed 's/\"//g' | sed 's/,//g' | awk '{print \$NF}')
    else
      tag="NA"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      nextclade: \$(nextclade --version | awk '{print \$NF}')
      tag: \$tag
      container: ${task.container}
    END_VERSIONS
  """
}
