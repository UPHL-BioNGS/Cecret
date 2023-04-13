process download {
  tag        "${sra}"
  publishDir "${params.outdir}", mode: 'copy'
  container  'quay.io/biocontainers/pandas:1.1.5'

  //#UPHLICA maxForks 10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '45m'

  input:
  val(sra)

  output:
  tuple val(sra), file("sra/paired/${sra}*fastq.gz"), val("paired"), optional: true, emit: paired
  tuple val(sra), file("sra/single/${sra}*fastq.gz"), val("single"), optional: true, emit: single

  shell:
  '''
    mkdir -p sra/{paired,single}

    sra=!{sra}

    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${sra:0:6}/0${sra: -2}/!{sra}/!{sra}_1.fastq.gz || \
        wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${sra:0:6}/0${sra: -2}/!{sra}/!{sra}.fastq.gz

    if [ -f !{sra}_1.fastq.gz ]
    then
        wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${sra:0:6}/0${sra: -2 }/!{sra}/!{sra}_2.fastq.gz
        mv !{sra}_*.fastq.gz sra/paired/.
    elif [ -f !{sra}.fastq.gz ]
    then
        mv !{sra}.fastq.gz sra/single/.
    else
        echo "Could not download file for SRA accession !{sra}"
    fi
  '''
}