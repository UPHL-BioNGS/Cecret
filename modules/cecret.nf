//# some fastas are created with the header of >reference, so this changes the header
process fasta_prep {
  tag        "${fasta}"
  //# nothing to publish in publishDir
  container  'quay.io/biocontainers/pandas:1.1.5'

  //#UPHLICA maxForks 10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '45m'

  when:
  fasta != null

  input:
  tuple val(sample), file(fasta)

  output:
  path "fasta_prep/${fasta}", optional: true, emit: fastas

  shell:
  '''
    mkdir -p fasta_prep

    echo ">!{sample}" > fasta_prep/!{fasta}
    grep -v ">" !{fasta} | fold -w 75 >> fasta_prep/!{fasta}
  '''
}

process unzip {
  tag        "${fasta}"
  //# nothing to publish in publishDir
  container  'quay.io/biocontainers/pandas:1.1.5'

  //#UPHLICA maxForks 10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '45m'

  when:
  params.nextclade || params.msa == 'nextalign'

  input:
  file(input)

  output:
  path "dataset", emit: dataset

  shell:
  '''
    mkdir dataset
    
    unzip !{input} -d dataset
  '''
}

process summary {
  tag        "${sample}"
  publishDir "${params.outdir}", mode: 'copy'
  container  'quay.io/biocontainers/pandas:1.1.5'

  //#UPHLICA maxForks 10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '45m'

  input:
  tuple val(sample), file(files), val(versions), path(multiqc), file(summary_script)

  output:
  path "summary/${sample}.summary.csv", emit: summary_file

  shell:
  '''
    mkdir -p summary

    exit 1

    python !{summary_script}
  '''
}
