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
  tag        "Creating summary files"
  publishDir "${params.outdir}", mode: 'copy'
  container  'quay.io/biocontainers/pandas:1.1.5'

  //#UPHLICA maxForks 10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '45m'

  input:
  tuple file(files), val(versions), file(multisample), path(multiqc), file(summary_script), file(fasta)

  output:
  path "cecret_results.{csv,txt}", emit: summary_file

  shell:
  '''
    echo "!{versions}" | cut -f 1,3,5,7,9,11  -d ',' | sed 's/\\[//g' | sed 's/\\]//g' | sed 's/, /,/g' >  versions.csv
    echo "!{versions}" | cut -f 2,4,6,8,10,12 -d ',' | sed 's/\\[//g' | sed 's/\\]//g' | sed 's/, /,/g' | awk '{($1=$1); print $0}' >> versions.csv

    grep -h ^FREADS *_ampliconstats.txt > ampliconstats.summary || echo "No ampliconstats file"

    if [ -s "vadr.vadr.sqa" ] ; then tail -n +2 "vadr.vadr.sqa" | grep -v "#-" | tr -s '[:blank:]' ',' > vadr.csv ; fi

    python !{summary_script}
  '''
}
