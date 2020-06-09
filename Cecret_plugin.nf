#!/usr/bin/env nextflow

println("Currently using the Cecret workflow adjusted for use after UPHL artic-Illumina hybrid library prep and pipeline")
println("v.20200912")

//# nextflow run /home/eriny/sandbox/Cecret/Cecret_plugin.nf -c /home/eriny/sandbox/Cecret/singularity.nextflow.config
//# To be used with the ivar container staphb/ivar:1.2.2_artic20200528, this includes all artic and reference files, plus the index files are pre-indexed
//# emacs covid_samples.txt where accession\tsubmission\tcollection_date
//# ~/sandbox/UPHL/COVID/files_for_submission.sh $(pwd)

maxcpus = Runtime.runtime.availableProcessors()
submission_script=workflow.projectDir + "/files_for_submission.sh"

params.artic_version = 'V3'

params.primer_bed = file("/artic-ncov2019/primer_schemes/nCoV-2019/${params.artic_version}/nCoV-2019.bed")
params.reference_genome = file("/home/eriny/src/artic-ncov2019/primer_schemes/nCoV-2019/${params.artic_version}/nCoV-2019.reference.fasta")
params.gff_file = file("/reference/GCF_009858895.2_ASM985889v3_genomic.gff")
params.amplicon_bed = file("/artic-ncov2019/primer_schemes/nCoV-2019/${params.artic_version}/nCoV-2019_amplicon.bed")

params.outdir = workflow.launchDir
params.log_directory = params.outdir + '/logs'
println("The files and directory for results is " + params.outdir)

params.sample_file = file(params.outdir + '/covid_samples.txt' )
if (!params.sample_file.exists()) exit 1, println "FATAL: ${params.sample_file} could not be found!\nPlease include a file name covid_samples.txt with the sample_id\tsubmission_id\tcollection_date at ${params.outdir}"
  else println "List of COVID19 samples: " + params.sample_file

samples = params.sample_file.readLines()[0]
samples_join = samples.join('|')

Channel
  .fromFilePairs("${params.outdir}/Sequencing_reads/QCed/*_clean_PE{1,2}.fastq", size: 2 )
  .ifEmpty{ exit 1, println("No cleaned COVID reads were found at ${params.outdir}/Sequencing_reads/QCed") }
  .set { clean_reads }

process bwa {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus maxcpus

  beforeScript 'mkdir -p covid/bwa logs/bwa_covid'

  input:
  set val(sample), file(reads) from clean_reads

  output:
  tuple sample, file("covid/bwa/${sample}.sorted.bam") into bams, bams2, bams3, bams4, bams5, bams6
  file("covid/bwa/${sample}.sorted.bam.bai") into bais
  file("logs/bwa_covid/${sample}.${workflow.sessionId}.log")
  file("logs/bwa_covid/${sample}.${workflow.sessionId}.err")

  when:
  sample =~ samples_join

  shell:
  '''
    log_file=logs/bwa_covid/!{sample}.!{workflow.sessionId}.log
    err_file=logs/bwa_covid/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "bwa $(bwa 2>&1 | grep Version )" >> $log_file
    samtools --version >> $log_file

    # bwa mem command
    bwa mem -t !{maxcpus} !{params.reference_genome} !{reads[0]} !{reads[1]} 2>> $err_file | \
      samtools sort 2>> $err_file | \
      samtools view -F 4 -o covid/bwa/!{sample}.sorted.bam 2>> $err_file >> $log_file

    # indexing the bams
    samtools index covid/bwa/!{sample}.sorted.bam 2>> $err_file >> $log_file
  '''
 }

process ivar_trim {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus 1

  beforeScript 'mkdir -p covid/trimmed logs/ivar_trim'

  input:
  set val(sample), file(bam) from bams

  output:
  tuple sample, file("covid/trimmed/${sample}.primertrim.bam") into trimmed_bams
  file("logs/ivar_trim/${sample}.${workflow.sessionId}.log")
  file("logs/ivar_trim/${sample}.${workflow.sessionId}.err")

  shell:
  '''
    log_file=logs/ivar_trim/!{sample}.!{workflow.sessionId}.log
    err_file=logs/ivar_trim/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    ivar version >> $log_file

    # trimming the reads
    ivar trim -e -i !{bam} -b !{params.primer_bed} -p covid/trimmed/!{sample}.primertrim 2>> $err_file >> $log_file
  '''
}

process samtools_sort {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus 1

  beforeScript 'mkdir -p covid/sorted logs/samtools_sort'

  input:
  set val(sample), file(bam) from trimmed_bams

  output:
  tuple sample, file("covid/sorted/${sample}.primertrim.sorted.bam") into sorted_bams, sorted_bams2, sorted_bams3, sorted_bams4, sorted_bams5
  file("covid/sorted/${sample}.primertrim.sorted.bam.bai") into sorted_bais
  file("logs/samtools_sort/${sample}.${workflow.sessionId}.log")
  file("logs/samtools_sort/${sample}.${workflow.sessionId}.err")

  shell:
  '''
    log_file=logs/samtools_sort/!{sample}.!{workflow.sessionId}.log
    err_file=logs/samtools_sort/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    # sorting and indexing the trimmed bams
    samtools sort !{bam} -o covid/sorted/!{sample}.primertrim.sorted.bam 2>> $err_file >> $log_file
    samtools index covid/sorted/!{sample}.primertrim.sorted.bam 2>> $err_file >> $log_file
  '''
}

process ivar_variants {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus 1

  beforeScript 'mkdir -p covid/variants logs/ivar_variants'

  input:
  set val(sample), file(bam) from sorted_bams

  output:
  tuple sample, file("covid/variants/${sample}.variants.tsv") into variants
  file("logs/ivar_variants/${sample}.${workflow.sessionId}.log")
  file("logs/ivar_variants/${sample}.${workflow.sessionId}.err")

  shell:
  '''
    log_file=logs/ivar_variants/!{sample}.!{workflow.sessionId}.log
    err_file=logs/ivar_variants/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file
    ivar version >> $log_file

    samtools mpileup -A -d 600000 -B -Q 0 --reference !{params.reference_genome} !{bam} 2>> $err_file | \
      ivar variants -p covid/variants/!{sample}.variants -q 20 -t 0.6 -r !{params.reference_genome} -g !{params.gff_file} 2>> $err_file >> $log_file
  '''
}

process ivar_consensus {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus 1

  beforeScript 'mkdir -p covid/consensus logs/ivar_consensus'

  input:
  set val(sample), file(bam) from sorted_bams2

  output:
  tuple sample, file("covid/consensus/${sample}.consensus.fa") into consensus, consensus2, consensus3
  file("logs/ivar_consensus/${sample}.${workflow.sessionId}.log")
  file("logs/ivar_consensus/${sample}.${workflow.sessionId}.err")

  shell:
  '''
    log_file=logs/ivar_consensus/!{sample}.!{workflow.sessionId}.log
    err_file=logs/ivar_consensus/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file
    ivar version >> $log_file


    samtools mpileup -A -d 6000000 -B -Q 0 --reference !{params.reference_genome} !{bam} 2>> $err_file | \
      ivar consensus -t 0.6 -p covid/consensus/!{sample}.consensus -n N 2>> $err_file >> $log_file
  '''
}

consensus
  .combine(bams2, by: 0)
  .set { combined_bams_consensus }

process quast {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus 1
  errorStrategy 'ignore'

  beforeScript 'mkdir -p covid/quast logs/quast_covid'

  input:
  set val(sample), file(fasta), file(bam) from combined_bams_consensus

  output:
  tuple sample, file("covid/quast/${sample}.report.txt") into quast_results
  path("covid/quast/${sample}")
  file("logs/quast_covid/${sample}.${workflow.sessionId}.log")
  file("logs/quast_covid/${sample}.${workflow.sessionId}.err")

  shell:
    '''
      log_file=logs/quast_covid/!{sample}.!{workflow.sessionId}.log
      err_file=logs/quast_covid/!{sample}.!{workflow.sessionId}.err

      date | tee -a $log_file $err_file > /dev/null
      quast.py --version >> $log_file


      quast.py !{fasta} -r !{params.reference_genome} --ref-bam !{bam} --output-dir covid/quast/!{sample} 2>> $err_file >> $log_file
      cp covid/quast/!{sample}/report.txt covid/quast/!{sample}.report.txt
    '''
}

bams3
  .combine(sorted_bams3, by: 0)
  .set { combined_bams3 }

process samtools_stats {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus 1

  beforeScript 'mkdir -p covid/samtools_stats/bwa covid/samtools_stats/sort logs/samtools_stats'

  input:
  set val(sample), file(bam), file(sorted_bam) from combined_bams3

  output:
  file("covid/samtools_stats/bwa/${sample}.stats.txt") into samtools_stats_results
  file("covid/samtools_stats/sort/${sample}.stats.trim.txt")
  file("logs/samtools_stats/${sample}.${workflow.sessionId}.log")
  file("logs/samtools_stats/${sample}.${workflow.sessionId}.err")

  shell:
  '''
    log_file=logs/samtools_stats/!{sample}.!{workflow.sessionId}.log
    err_file=logs/samtools_stats/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file


    samtools stats !{bam} > covid/samtools_stats/bwa/!{sample}.stats.txt 2>> $err_file
    samtools stats !{sorted_bam} > covid/samtools_stats/sort/!{sample}.stats.trim.txt 2>> $err_file
  '''
}

bams4
  .combine(sorted_bams4, by: 0)
  .set { combined_bams4 }

process samtools_coverage {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus 1

  beforeScript 'mkdir -p covid/samtools_coverage/bwa covid/samtools_coverage/sort logs/samtools_coverage'

  input:
  set val(sample), file(bwa), file(sorted) from combined_bams4

  output:
  file("covid/samtools_coverage/bwa/${sample}.cov.txt") into samtools_coverage_results
  file("covid/samtools_coverage/bwa/${sample}.cov.hist")
  file("covid/samtools_coverage/sort/${sample}.cov.trim.txt")
  file("covid/samtools_coverage/sort/${sample}.cov.trim.hist")
  file("logs/samtools_coverage/${sample}.${workflow.sessionId}.log")
  file("logs/samtools_coverage/${sample}.${workflow.sessionId}.err")

  shell:
  '''
    log_file=logs/samtools_coverage/!{sample}.!{workflow.sessionId}.log
    err_file=logs/samtools_coverage/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file


    samtools coverage !{bwa} -m -o covid/samtools_coverage/bwa/!{sample}.cov.hist 2>> $err_file >> $log_file
    samtools coverage !{bwa} -o covid/samtools_coverage/bwa/!{sample}.cov.txt 2>> $err_file >> $log_file
    samtools coverage !{sorted} -m -o covid/samtools_coverage/sort/!{sample}.cov.trim.hist 2>> $err_file >> $log_file
    samtools coverage !{sorted} -o covid/samtools_coverage/sort/!{sample}.cov.trim.txt 2>> $err_file >> $log_file
  '''
}

process bedtools {
  publishDir "${params.outdir}", mode: 'copy'
  tag "bedtools"
  echo true
  cpus 1

  beforeScript 'mkdir -p covid/bedtools logs/bedtools'

  input:
  file(bwa) from bams5.collect()
  file(sort) from sorted_bams5.collect()
  file(bai) from bais.collect()
  file(sorted_bai) from sorted_bais.collect()

  output:
  file("covid/bedtools/multicov.txt") into bedtools_results
  file("logs/bedtools/multicov.${workflow.sessionId}.log")
  file("logs/bedtools/multicov.${workflow.sessionId}.err")

  shell:
  '''
    log_file=logs/bedtools/multicov.!{workflow.sessionId}.log
    err_file=logs/bedtools/multicov.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    bedtools --version >> $log_file


    echo "primer" $(ls *bam) | tr ' ' '\t' > covid/bedtools/multicov.txt
    bedtools multicov -bams $(ls *bam) -bed !{params.amplicon_bed} | cut -f 4,6- 2>> $err_file >> covid/bedtools/multicov.txt
  '''
}

process summary {
  publishDir "${params.outdir}", mode: 'copy', overwrite: true
  tag "summary"
  echo true
  cpus 1

  beforeScript 'mkdir -p covid logs/summary'

  input:
  file(bwa) from bams6.collect()
  file(consensus) from consensus2.collect()
  file(quast) from quast_results.collect()
  file(stats) from samtools_stats_results.collect()
  file(coverage) from samtools_coverage_results.collect()
  file(multicov) from bedtools_results.collect()

  output:
  file("covid/summary.txt") into final_summary
  file("logs/summary/summary.${workflow.sessionId}.log")
  file("logs/summary/summary.${workflow.sessionId}.err")

  shell:
  '''
    log_file=logs/summary/summary.!{workflow.sessionId}.log
    err_file=logs/summary/summary.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null


    echo "sample,%_Human_reads,degenerate_check,coverage,depth,failed_amplicons,num_N" > covid/summary.txt

    while read line
    do
      sample=$(grep $line !{params.outdir}/run_results.txt | cut -f 2 | head -n 1 )
      find_test=$(find !{params.outdir}/Sequencing_reads/QCed/. -iname "$line*" | head -n 1 )
      if [ -n "$find_test" ]
      then
        human_reads=$(grep "Homo" !{params.outdir}/blobtools/$line*100.blobplot.stats.txt | cut -f 13 ) 2>> $err_file
        if [ -z "$human_reads" ] ; then human_reads="none" ; fi

        degenerate=$(grep -f ~/degenerate.txt $line*.consensus.fa | grep -v ">" | wc -l ) 2>> $err_file
        if [ -z "$degenerate" ] ; then degenerate="none" ; fi

        cov_and_depth=($(cut -f 6,7 $line*.cov.txt | tail -n 1)) 2>> $err_file
        if [ -z "${cov_and_depth[0]}" ] ; then cov_and_depth=(0 0) ; fi

        bedtools_column=$(head -n 1 multicov.txt | tr '\t' '\n' | grep -n $line | grep -v primertrim | cut -f 1 -d ":" | head -n 1 ) 2>> $err_file
        amp_fail=$(cut -f $bedtools_column multicov.txt | awk '{{ if ( $1 < 20 ) print $0 }}' | wc -l ) 2>> $err_file
        if [ -z "$amp_fail" ] ; then amp_fail=0 ; fi

        num_of_N=$(grep -o 'N' $line*.consensus.fa | wc -l ) 2>> $err_file
        if [ -z "$num_of_N" ] ; then num_of_N=0 ; fi

        echo "$sample,$human_reads,$degenerate,${cov_and_depth[0]},${cov_and_depth[1]},$amp_fail,$num_of_N" >> covid/summary.txt
      fi

    done < <(cat !{params.sample_file} | awk '{print $1}')
  '''
}

process file_submission {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus 1

  beforeScript 'mkdir -p covid/files_for_submission logs/submission'

  input:
  file(summary) from consensus3.collect()

  output:

  shell:
  '''
  log_file=logs/submission/submission.!{workflow.sessionId}.log
  err_file=logs/submission/submission.!{workflow.sessionId}.err

  !{submission_script} !{params.outdir}
  '''
}

process multiqc {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus 1

  beforeScript 'mkdir -p covid/multiqc logs/multiqc'

  input:
  file(file) from final_summary.collect()

  output:
    path("covid/multiqc/multiqc_data")
    file("covid/multiqc/multiqc_report.html")

  shell:
    '''
    log_file=logs/multiqc/multiqc.!{workflow.sessionId}.log
    err_file=logs/multiqc/multiqc.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    multiqc --version >> $log_file

    multiqc -f \
        --outdir covid/multiqc \
        !{params.outdir}/covid
        2>> $err_file | tee -a $log_file
    '''
}

workflow.onComplete {
    println("Pipeline completed at: $workflow.complete")
    println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
}
