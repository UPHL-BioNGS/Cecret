//# some fastas are created with the header of >reference, so this changes the header
process fasta_prep {
  tag "${fasta}"

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

process summary {
  tag "${sample}"

  input:
  tuple val(sample), file(files), file(summary_script)

  output:
  path "summary/${sample}.summary.csv", emit: summary_file

  shell:
  '''
    mkdir -p summary

    python !{summary_script}    
  '''
}
