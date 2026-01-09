process DATASETS {
  tag           "${accession}"
  // because there's no way to specify threads
  label         "process_low"
  container     'staphb/ncbi-datasets:18.13.0'

  input:
  val(accession)

  output:
  path "genomes/*fasta", emit: fasta, optional: true
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args   ?: "virus"
  """
    mkdir -p genomes

    datasets \
      download \
      ${args} \
      genome \
      accession \
      ${accession} \
      --filename ncbi_dataset.zip

    unzip ncbi_dataset.zip

    cp ncbi_dataset/data/genomic.fna genomes/${accession}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      datasets: \$(datasets --version | awk '{print \$NF}')
    END_VERSIONS
  """
}
