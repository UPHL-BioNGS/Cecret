//# some fastas are created with the header of >reference, so this changes the header
process PREP {
  tag        "${fasta}"
  container  'staphb/pandas:2.3.3'
  label      "process_single"

  input:
  tuple val(meta), file(fasta)

  output:
  path "fasta_prep/*", optional: true, emit: fastas
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p fasta_prep

  echo ">${prefix}" > fasta_prep/${prefix}.fasta
  grep -v ">" ${fasta} | fold -w 75 >> fasta_prep/${prefix}.fasta

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    prep: NA
    container: ${task.container}
  END_VERSIONS
  """
}

process SUMMARY {
  tag        "Creating summary files"
  label      "process_low"
  container  'staphb/pandas:2.3.3'

  input:
  tuple file(files), file(script), file(multiqc)

  output:
  path "cecret_results.{csv,txt}", emit: summary_file

  when:
  task.ext.when == null || task.ext.when

  script:
  def multiqc_files = multiqc.join(" ")
  """
    mkdir multiqc_data
    for file in ${multiqc_files}
    do
      if [ -f "\$file" ]; then mv \$file multiqc_data/. ; fi
    done

    if [ -n "\$(find . -iname *ampliconstats.txt | head -n 1)" ] 
    then
      cat *_ampliconstats.txt | grep -h ^FREADS > ampliconstats.summary
    else
      touch ampliconstats.summary
    fi

    if [ -s "vadr.vadr.sqa" ] ; then tail -n +2 "vadr.vadr.sqa" | grep -v "#-" | tr -s '[:blank:]' ',' > vadr.csv ; fi

    python3 combine_results.py ${params.minimum_depth} ${workflow.manifest.version}
  """
}

process UNZIP {
  tag        "unzipping nextclade dataset"
  label      "process_single"
  container  'staphb/ncbi-datasets:18.9.0'

  input:
  file(input)

  output:
  path "dataset", emit: dataset
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  """
  mkdir dataset
  unzip ${input} -d dataset

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    unzip: NA
    container: ${task.container}
  END_VERSIONS
  """
}
