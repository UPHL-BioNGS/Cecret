#!/usr/bin/env nextflow

println("UPHL Reference Free Workflow v.20200605")

//# nextflow run ~/sandbox/UPHL/UPHL_reference_free_docker.nf
// nextflow run ~/sandbox/UPHL/UPHL_reference_free_docker.nf -c ~/sandbox/UPHL/URF_scripts/singularity.nextflow.config

params.outdir = workflow.launchDir
params.blast_db = "/home/Bioinformatics/Data/blastdb/"
params.blast_db_container = "/blast/blastdb"


maxcpus = Runtime.runtime.availableProcessors()
println("The maximum number of CPUS used in this workflow is ${maxcpus}")
if ( maxcpus < 5 ) {
  medcpus = maxcpus
} else {
  medcpus = 5
}

maxmem = Math.round(Runtime.runtime.totalMemory() / 10241024)
println("The maximum amount of memory used in this workflow is ${maxmem}")

Channel
  .fromFilePairs(["${params.outdir}/Sequencing_reads/Raw/*_R{1,2}_001.fastq.gz", "${params.outdir}/Sequencing_reads/Raw/*_{1,2}.fastq" ], size: 2 )
  .map{ reads -> [reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1]] }
  .ifEmpty{ exit 1, println("No paired fastq or fastq.gz files were found at ${params.outdir}/Sequencing_reads/Raw") }
  .into { fastq_reads; fastq_reads2 }

process seqyclean {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo true
  cpus 1

  beforeScript 'mkdir -p Sequencing_reads/QCed logs/seqyclean'

  input:
  set val(sample), file(reads) from fastq_reads

  output:
  tuple sample, file("Sequencing_reads/QCed/${sample}_clean_PE{1,2}.fastq") into clean_reads, clean_reads2
  file("Sequencing_reads/QCed/${sample}_clean_SE.fastq")
  file("Sequencing_reads/QCed/${sample}_clean_SummaryStatistics.{txt,tsv}")
  file("logs/seqyclean/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/seqyclean/!{sample}.!{workflow.sessionId}.log
    err_file=logs/seqyclean/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "seqyclean version: $(seqyclean -h | grep Version)" >> $log_file

    seqyclean -minlen 25 -qual -c !{params.seqyclean_contaminant_file} -1 !{reads[0]} -2 !{reads[1]} -o Sequencing_reads/QCed/!{sample}_clean 2>> $err_file >> $log_file
  '''
}

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

// fastq_reads2
//   .combine(clean_reads2, by: 0)
//   .into { raw_clean_reads; raw_clean_reads2 }
//
// process fastqc {
//   publishDir "${params.outdir}", mode: 'copy'
//   tag "$sample"
//   echo true
//   cpus 1
//
//   beforeScript 'mkdir -p fastqc logs/fastqc'
//
//   input:
//   set val(sample), file(raw), file(clean) from raw_clean_reads
//
//   output:
//   path("fastqc/")
//   tuple sample, env(raw_1), env(raw_2), env(clean_1), env(clean_2) into fastqc_results
//   file("logs/fastqc/${sample}.${workflow.sessionId}.{log,err}")
//
//   shell:
//   '''
//     log_file=logs/fastqc/!{sample}.!{workflow.sessionId}.log
//     err_file=logs/fastqc/!{sample}.!{workflow.sessionId}.err
//
//     # time stamp + capturing tool versions
//     date | tee -a $log_file $err_file > /dev/null
//     fastqc --version >> $log_file
//
//     fastqc --outdir fastqc --threads !{task.cpus} !{sample}*.fastq* 2>> $err_file >> $log_file
//
//     raw_1=$(unzip -p fastqc/!{raw[0].simpleName}*fastqc.zip */fastqc_data.txt | grep "Total Sequences" | awk '{ print $3 }' )
//     raw_2=$(unzip -p fastqc/!{raw[1].simpleName}*fastqc.zip */fastqc_data.txt | grep "Total Sequences" | awk '{ print $3 }' )
//     clean_1=$(unzip -p fastqc/!{clean[0].simpleName}*fastqc.zip */fastqc_data.txt | grep "Total Sequences" | awk '{ print $3 }' )
//     clean_2=$(unzip -p fastqc/!{clean[1].simpleName}*fastqc.zip */fastqc_data.txt | grep "Total Sequences" | awk '{ print $3 }' )
//   '''
// }
//
//
//
// contigs
//   .combine(mash_organism, by: 0)
//   .set { contig_organism }
//
//
//
// process quast {
//   publishDir "${params.outdir}", mode: 'copy'
//   tag "$sample"
//   echo true
//   cpus medcpus
//
//   beforeScript 'mkdir -p quast logs/quast'
//
//   input:
//   set val(sample), file(contig) from contigs4
//
//   output:
//   file("quast/${sample}/icarus.html")
//   path("quast/${sample}/{basic_stats,icarus_viewers}")
//   file("quast/${sample}/report.{html,pdf,tex,txt,tsv}")
//   file("quast/${sample}/transposed_report.{tex,tsv,txt}")
//   file("logs/quast/${sample}.${workflow.sessionId}.{log,err}")
//
//   shell:
//   '''
//     log_file=logs/quast/!{sample}.!{workflow.sessionId}.log
//     err_file=logs/quast/!{sample}.!{workflow.sessionId}.err
//
//     # time stamp + capturing tool versions
//     date | tee -a $log_file $err_file > /dev/null
//     quast.py --version >> $log_file
//
//     quast.py !{contig} --output-dir quast/!{sample} --threads !{task.cpus} 2>> $err_file | tee -a $log_file
//   '''
// }
//
//
//
// process blastn {
//   publishDir "${params.outdir}", mode: 'copy'
//   tag "$sample"
//   echo true
//   cpus medcpus
//
//   beforeScript 'mkdir -p blast logs/blastn'
//
//   input:
//   set val(sample), file(contigs) from contigs7
//
//   output:
//   tuple sample, file("blast/${sample}.tsv") into blast
//   file("logs/blastn/${sample}.${workflow.sessionId}.{log,err}")
//
//   shell:
//   '''
//     log_file=logs/blastn/!{sample}.!{workflow.sessionId}.log
//     err_file=logs/blastn/!{sample}.!{workflow.sessionId}.err
//
//     # time stamp + capturing tool versions
//     date | tee -a $log_file $err_file > /dev/null
//     blastn -version >> $log_file
//     echo "The blastdb location is $BLASTDB" >> $log_file
//
//     blastn -query !{contigs} \
//       -out blast/!{sample}.tsv \
//       -num_threads !{task.cpus} \
//       -db !{params.blast_db_container}/nt \
//       -outfmt '6 qseqid staxids bitscore std' \
//       -max_target_seqs 10 \
//       -max_hsps 1 \
//       -evalue 1e-25 2>> $err_file >> $log_file
//    '''
// }
//
// index
//   .combine(clean_reads5, by:0)
//   .combine(contigs8, by:0)
//   .set { index_clean_contigs }
//
// process bwa {
//   publishDir "${params.outdir}", mode: 'copy'
//   tag "$sample"
//   echo true
//   cpus maxcpus
//
//   beforeScript 'mkdir -p bwa logs/bwa'
//
//   input:
//   set val(sample), file(index), file(reads), file(contigs) from index_clean_contigs
//
//   output:
//   tuple sample, file("bwa/${sample}.sorted.bam"), file("bwa/${sample}.sorted.bam.bai") into bams
//   file("logs/bwa/${sample}.${workflow.sessionId}.{log,err}")
//
//   shell:
//   '''
//     log_file=logs/bwa/!{sample}.!{workflow.sessionId}.log
//     err_file=logs/bwa/!{sample}.!{workflow.sessionId}.err
//
//     # time stamp + capturing tool versions
//     date | tee -a $log_file $err_file > /dev/null
//     echo "bwa $(bwa 2>&1 | grep Version )" >> $log_file
//     samtools --version >> {output.log}.log
//
//     bwa mem -t !{task.cpus} !{contigs} !{reads[0]} !{reads[1]} 2>> $err_file | \
//       samtools sort -o bwa/!{sample}.sorted.bam 2>> $err_file >> $log_file
//
//     samtools index bwa/!{sample}.sorted.bam 2>> $err_file >> $log_file
//
//     exit 1
//   '''
// }
//
// blast
//   .combine(bams, by:0)
//   .combine(contigs2, by:0)
//   .set { blastn_bwa_contigs }
//
// process blobtools_create {
//   publishDir "${params.outdir}", mode: 'copy'
//   tag "$sample"
//   echo true
//   cpus 1
//
//   beforeScript 'mkdir -p blobtools logs/blobtools_create'
//
//   input:
//   set val(sample), file(blast), file(bwa), file(contigs) from blastn_bwa_contigs
//
//   output:
//   file("blobtools/${sample}.${sample}.sorted.bam.cov")
//   tuple sample, file("blobtools/${sample}.blobDB.json") into create, create2
//   file("logs/blobtools_create/${sample}.${workflow.sessionId}.{log,err}")
//
//   shell:
//   '''
//     log_file=logs/blobtools_create/!{sample}.!{workflow.sessionId}.log
//     err_file=logs/blobtools_create/!{sample}.!{workflow.sessionId}.err
//
//     # time stamp + capturing tool versions
//     date | tee -a $log_file $err_file > /dev/null
//     echo "blobtools version $(blobtools -v)" >> $log_file
//
//     blobtools create -o blobtools/!{sample} -i !{contigs} -b !{bam} -t !{blast} 2>> $err_file >> $log_file
//
//     exit 1
//   '''
// }
//
// process blobtools_view {
//   publishDir "${params.outdir}", mode: 'copy'
//   tag "$sample"
//   echo true
//   cpus 1
//
//   beforeScript 'mkdir -p logs/blobtools_view'
//
//   input:
//   set val(sample), file(json) from create
//
//   output:
//   tuple sample, file("blobtools/${sample}.blobDB.table.txt") into view
//   file("logs/blobtools_view/${sample}.${workflow.sessionId}.{log,err}")
//
//   shell:
//   '''
//     log_file=logs/blobtools_view/!{sample}.!{workflow.sessionId}.log
//     err_file=logs/blobtools_view/!{sample}.!{workflow.sessionId}.err
//
//     # time stamp + capturing tool versions
//     date | tee -a $log_file $err_file > /dev/null
//     echo "blobtools version $(blobtools -v)" >> $log_file
//
//     blobtools view -i !{json} -o blobtools/ 2>> $err_file >> $log_file
//
//     exit 1
//   '''
// }
//
// create2
//   .combine(view, by:0)
//   .set { create_view }
//
// process blobtools_plot {
//   publishDir "${params.outdir}", mode: 'copy'
//   tag "$sample"
//   echo true
//   cpus 1
//
//   beforeScript 'mkdir -p logs/blobtools_plot'
//
//   input:
//   set val(sample), file(json), file(view) from create_view
//
//   output:
//   file("blobtools/${sample}.blobDB.json.bestsum.species.p8.span.100.blobplot.bam0.png")
//   file("blobtools/${sample}.blobDB.json.bestsum.species.p8.span.100.blobplot.read_cov.bam0.png")
//   file("blobtools/${sample}.blobDB.json.bestsum.species.p8.span.100.blobplot.stats.txt")
//   tuple sample, env(result) into blobtools_results
//   file("logs/blobtools_plot/${sample}.${workflow.sessionId}.{log,err}")
//
//   shell:
//   '''
//     log_file=logs/blobtools_plot/!{sample}.!{workflow.sessionId}.log
//     err_file=logs/blobtools_plot/!{sample}.!{workflow.sessionId}.err
//
//     # time stamp + capturing tool versions
//     date | tee -a $log_file $err_file > /dev/null
//     echo "blobtools version $(blobtools -v)" >> $log_file
//
//     blobtools plot -i !{cov} -o blobtools/ -r species --format png 2>> $err_file >> $log_file
//
//     result=$(grep -v ^"#" $out/blobtools/$sample.blobDB.json.bestsum.species.p8.span.100.blobplot.stats.txt | grep -v ^"all" | head -n 1 | tr ' ' '_' | cut -f 1,13)
//
//     exit 1
//   '''
// }
//
// //process pangolin
//
// fastqc_results
//   .combine(mash_results, by:0)
//   .combine(cg_pipeline_results, by:0)
//   .combine(seqsero2_results, by:0)
//   .combine(abricate_results, by:0)
//   .combine(blobtools_results, by:0)
//   .combine(mlst_results, by:0)
//   .set { workflow_results }
//
// workflow_results.view()
//
// // process multiqc {
// //   publishDir "${params.outdir}", mode: 'copy'
// //   tag "$sample"
// //   echo true
// //   cpus 1
// //
// //   beforeScript 'mkdir -p multiqc logs/multiqc'
// //
// //   input:
// //
// //   output:
// //   tuple sample, file("mash/${sample}_mashdist.txt")
// //   file("logs/multiqc/${sample}.${workflow.sessionId}.{log,err}")
// //
// //   shell:
// //   '''
// //     log_file=logs/multiqc/multiqc.!{workflow.sessionId}.log
// //     err_file=logs/multiqc/multiqc.!{workflow.sessionId}.err
// //
// //     # time stamp + capturing tool versions
// //     date | tee -a $log_file $err_file > /dev/null
// //     multiqc --version >> $log_file
// //
// //     multiqc -f \
// //       --outdir ./logs \
// //       --cl_config "prokka_fn_snames: True"  \
// //       !{params.outdir}/abricate_results/summary \
// //       !{params.outdir}/blobtools \
// //       !{params.outdir}/cg-pipeline \
// //       !{params.outdir}/fastqc \
// //       !{params.outdir}/mash \
// //       !{params.outdir}/prokka \
// //       !{params.outdir}/quast \
// //       !{params.outdir}/seqsero2 \
// //       !{params.outdir}/Sequencing_reads/QCed/*tsv \
// //       2>> $err_file | tee -a $log_file
// //
// //     exit 1
// //   '''
// // }
//
// process ivar_trim {
//   publishDir "${params.outdir}", mode: 'copy'
//   tag "$sample"
//   echo true
//   cpus 1
//
//   beforeScript 'mkdir -p covid/trimmed logs/ivar_trim'
//
//   input:
//   set val(sample), file(bam) from bams
//
//   output:
//   tuple sample, file("covid/trimmed/${sample}.primertrim.bam") into trimmed_bams
//   file("logs/ivar_trim/${sample}.${workflow.sessionId}.log")
//   file("logs/ivar_trim/${sample}.${workflow.sessionId}.err")
//
//   shell:
//   '''
//     log_file=logs/ivar_trim/!{sample}.!{workflow.sessionId}.log
//     err_file=logs/ivar_trim/!{sample}.!{workflow.sessionId}.err
//
//     # time stamp + capturing tool versions
//     date | tee -a $log_file $err_file > /dev/null
//     ivar version >> $log_file
//
//     # trimming the reads
//     ivar trim -e -i !{bam} -b !{params.primer_bed} -p covid/trimmed/!{sample}.primertrim 2>> $err_file >> $log_file
//   '''
// }
//
// process samtools_sort {
//   publishDir "${params.outdir}", mode: 'copy'
//   tag "$sample"
//   echo true
//   cpus 1
//
//   beforeScript 'mkdir -p covid/sorted logs/samtools_sort'
//
//   input:
//   set val(sample), file(bam) from trimmed_bams
//
//   output:
//   tuple sample, file("covid/sorted/${sample}.primertrim.sorted.bam") into sorted_bams, sorted_bams2, sorted_bams3, sorted_bams4, sorted_bams5
//   file("covid/sorted/${sample}.primertrim.sorted.bam.bai") into sorted_bais
//   file("logs/samtools_sort/${sample}.${workflow.sessionId}.log")
//   file("logs/samtools_sort/${sample}.${workflow.sessionId}.err")
//
//   shell:
//   '''
//     log_file=logs/samtools_sort/!{sample}.!{workflow.sessionId}.log
//     err_file=logs/samtools_sort/!{sample}.!{workflow.sessionId}.err
//
//     # time stamp + capturing tool versions
//     date | tee -a $log_file $err_file > /dev/null
//     samtools --version >> $log_file
//
//     # sorting and indexing the trimmed bams
//     samtools sort !{bam} -o covid/sorted/!{sample}.primertrim.sorted.bam 2>> $err_file >> $log_file
//     samtools index covid/sorted/!{sample}.primertrim.sorted.bam 2>> $err_file >> $log_file
//   '''
// }
//
// process ivar_variants {
//   publishDir "${params.outdir}", mode: 'copy'
//   tag "$sample"
//   echo true
//   cpus 1
//
//   beforeScript 'mkdir -p covid/variants logs/ivar_variants'
//
//   input:
//   set val(sample), file(bam) from sorted_bams
//
//   output:
//   tuple sample, file("covid/variants/${sample}.variants.tsv") into variants
//   file("logs/ivar_variants/${sample}.${workflow.sessionId}.log")
//   file("logs/ivar_variants/${sample}.${workflow.sessionId}.err")
//
//   shell:
//   '''
//     log_file=logs/ivar_variants/!{sample}.!{workflow.sessionId}.log
//     err_file=logs/ivar_variants/!{sample}.!{workflow.sessionId}.err
//
//     # time stamp + capturing tool versions
//     date | tee -a $log_file $err_file > /dev/null
//     samtools --version >> $log_file
//     ivar version >> $log_file
//
//     samtools mpileup -A -d 600000 -B -Q 0 --reference !{params.reference_genome} !{bam} 2>> $err_file | \
//       ivar variants -p covid/variants/!{sample}.variants -q 20 -t 0.6 -r !{params.reference_genome} -g !{params.gff_file} 2>> $err_file >> $log_file
//   '''
// }
//
// process ivar_consensus {
//   publishDir "${params.outdir}", mode: 'copy'
//   tag "$sample"
//   echo true
//   cpus 1
//
//   beforeScript 'mkdir -p covid/consensus logs/ivar_consensus'
//
//   input:
//   set val(sample), file(bam) from sorted_bams2
//
//   output:
//   tuple sample, file("covid/consensus/${sample}.consensus.fa") into consensus, consensus2, consensus3
//   file("logs/ivar_consensus/${sample}.${workflow.sessionId}.log")
//   file("logs/ivar_consensus/${sample}.${workflow.sessionId}.err")
//
//   shell:
//   '''
//     log_file=logs/ivar_consensus/!{sample}.!{workflow.sessionId}.log
//     err_file=logs/ivar_consensus/!{sample}.!{workflow.sessionId}.err
//
//     date | tee -a $log_file $err_file > /dev/null
//     samtools --version >> $log_file
//     ivar version >> $log_file
//
//
//     samtools mpileup -A -d 6000000 -B -Q 0 --reference !{params.reference_genome} !{bam} 2>> $err_file | \
//       ivar consensus -t 0.6 -p covid/consensus/!{sample}.consensus -n N 2>> $err_file >> $log_file
//   '''
// }
//
// consensus
//   .combine(bams2, by: 0)
//   .set { combined_bams_consensus }
//
// process quast {
//   publishDir "${params.outdir}", mode: 'copy'
//   tag "$sample"
//   echo true
//   cpus 1
//   errorStrategy 'ignore'
//
//   beforeScript 'mkdir -p covid/quast logs/quast_covid'
//
//   input:
//   set val(sample), file(fasta), file(bam) from combined_bams_consensus
//
//   output:
//   tuple sample, file("covid/quast/${sample}.report.txt") into quast_results
//   path("covid/quast/${sample}")
//   file("logs/quast_covid/${sample}.${workflow.sessionId}.log")
//   file("logs/quast_covid/${sample}.${workflow.sessionId}.err")
//
//   shell:
//     '''
//       log_file=logs/quast_covid/!{sample}.!{workflow.sessionId}.log
//       err_file=logs/quast_covid/!{sample}.!{workflow.sessionId}.err
//
//       date | tee -a $log_file $err_file > /dev/null
//       quast.py --version >> $log_file
//
//
//       quast.py !{fasta} -r !{params.reference_genome} --ref-bam !{bam} --output-dir covid/quast/!{sample} 2>> $err_file >> $log_file
//       cp covid/quast/!{sample}/report.txt covid/quast/!{sample}.report.txt
//     '''
// }
//
// bams3
//   .combine(sorted_bams3, by: 0)
//   .set { combined_bams3 }
//
// process samtools_stats {
//   publishDir "${params.outdir}", mode: 'copy'
//   tag "$sample"
//   echo true
//   cpus 1
//
//   beforeScript 'mkdir -p covid/samtools_stats/bwa covid/samtools_stats/sort logs/samtools_stats'
//
//   input:
//   set val(sample), file(bam), file(sorted_bam) from combined_bams3
//
//   output:
//   file("covid/samtools_stats/bwa/${sample}.stats.txt") into samtools_stats_results
//   file("covid/samtools_stats/sort/${sample}.stats.trim.txt")
//   file("logs/samtools_stats/${sample}.${workflow.sessionId}.log")
//   file("logs/samtools_stats/${sample}.${workflow.sessionId}.err")
//
//   shell:
//   '''
//     log_file=logs/samtools_stats/!{sample}.!{workflow.sessionId}.log
//     err_file=logs/samtools_stats/!{sample}.!{workflow.sessionId}.err
//
//     date | tee -a $log_file $err_file > /dev/null
//     samtools --version >> $log_file
//
//
//     samtools stats !{bam} > covid/samtools_stats/bwa/!{sample}.stats.txt 2>> $err_file
//     samtools stats !{sorted_bam} > covid/samtools_stats/sort/!{sample}.stats.trim.txt 2>> $err_file
//   '''
// }
//
// bams4
//   .combine(sorted_bams4, by: 0)
//   .set { combined_bams4 }
//
// process samtools_coverage {
//   publishDir "${params.outdir}", mode: 'copy'
//   tag "$sample"
//   echo true
//   cpus 1
//
//   beforeScript 'mkdir -p covid/samtools_coverage/bwa covid/samtools_coverage/sort logs/samtools_coverage'
//
//   input:
//   set val(sample), file(bwa), file(sorted) from combined_bams4
//
//   output:
//   file("covid/samtools_coverage/bwa/${sample}.cov.txt") into samtools_coverage_results
//   file("covid/samtools_coverage/bwa/${sample}.cov.hist")
//   file("covid/samtools_coverage/sort/${sample}.cov.trim.txt")
//   file("covid/samtools_coverage/sort/${sample}.cov.trim.hist")
//   file("logs/samtools_coverage/${sample}.${workflow.sessionId}.log")
//   file("logs/samtools_coverage/${sample}.${workflow.sessionId}.err")
//
//   shell:
//   '''
//     log_file=logs/samtools_coverage/!{sample}.!{workflow.sessionId}.log
//     err_file=logs/samtools_coverage/!{sample}.!{workflow.sessionId}.err
//
//     date | tee -a $log_file $err_file > /dev/null
//     samtools --version >> $log_file
//
//
//     samtools coverage !{bwa} -m -o covid/samtools_coverage/bwa/!{sample}.cov.hist 2>> $err_file >> $log_file
//     samtools coverage !{bwa} -o covid/samtools_coverage/bwa/!{sample}.cov.txt 2>> $err_file >> $log_file
//     samtools coverage !{sorted} -m -o covid/samtools_coverage/sort/!{sample}.cov.trim.hist 2>> $err_file >> $log_file
//     samtools coverage !{sorted} -o covid/samtools_coverage/sort/!{sample}.cov.trim.txt 2>> $err_file >> $log_file
//   '''
// }
//
// process bedtools {
//   publishDir "${params.outdir}", mode: 'copy'
//   tag "bedtools"
//   echo true
//   cpus 1
//
//   beforeScript 'mkdir -p covid/bedtools logs/bedtools'
//
//   input:
//   file(bwa) from bams5.collect()
//   file(sort) from sorted_bams5.collect()
//   file(bai) from bais.collect()
//   file(sorted_bai) from sorted_bais.collect()
//
//   output:
//   file("covid/bedtools/multicov.txt") into bedtools_results
//   file("logs/bedtools/multicov.${workflow.sessionId}.log")
//   file("logs/bedtools/multicov.${workflow.sessionId}.err")
//
//   shell:
//   '''
//     log_file=logs/bedtools/multicov.!{workflow.sessionId}.log
//     err_file=logs/bedtools/multicov.!{workflow.sessionId}.err
//
//     date | tee -a $log_file $err_file > /dev/null
//     bedtools --version >> $log_file
//
//
//     echo "primer" $(ls *bam) | tr ' ' '\t' > covid/bedtools/multicov.txt
//     bedtools multicov -bams $(ls *bam) -bed !{params.amplicon_bed} | cut -f 4,6- 2>> $err_file >> covid/bedtools/multicov.txt
//   '''
// }
//
// process summary {
//   publishDir "${params.outdir}", mode: 'copy', overwrite: true
//   tag "summary"
//   echo true
//   cpus 1
//
//   beforeScript 'mkdir -p covid logs/summary'
//
//   input:
//   file(bwa) from bams6.collect()
//   file(consensus) from consensus2.collect()
//   file(quast) from quast_results.collect()
//   file(stats) from samtools_stats_results.collect()
//   file(coverage) from samtools_coverage_results.collect()
//   file(multicov) from bedtools_results.collect()
//
//   output:
//   file("covid/summary.txt") into final_summary
//   file("logs/summary/summary.${workflow.sessionId}.log")
//   file("logs/summary/summary.${workflow.sessionId}.err")
//
//   shell:
//   '''
//     log_file=logs/summary/summary.!{workflow.sessionId}.log
//     err_file=logs/summary/summary.!{workflow.sessionId}.err
//
//     date | tee -a $log_file $err_file > /dev/null
//
//
//     echo "sample,%_Human_reads,degenerate_check,coverage,depth,failed_amplicons,num_N" > covid/summary.txt
//
//     while read line
//     do
//       sample=$(grep $line !{params.outdir}/run_results.txt | cut -f 2 | head -n 1 )
//       find_test=$(find !{params.outdir}/Sequencing_reads/QCed/. -iname "$line*" | head -n 1 )
//       if [ -n "$find_test" ]
//       then
//         human_reads=$(grep "Homo" !{params.outdir}/blobtools/$line*100.blobplot.stats.txt | cut -f 13 ) 2>> $err_file
//         if [ -z "$human_reads" ] ; then human_reads="none" ; fi
//
//         degenerate=$(grep -f ~/degenerate.txt $line*.consensus.fa | grep -v ">" | wc -l ) 2>> $err_file
//         if [ -z "$degenerate" ] ; then degenerate="none" ; fi
//
//         cov_and_depth=($(cut -f 6,7 $line*.cov.txt | tail -n 1)) 2>> $err_file
//         if [ -z "${cov_and_depth[0]}" ] ; then cov_and_depth=(0 0) ; fi
//
//         bedtools_column=$(head -n 1 multicov.txt | tr '\t' '\n' | grep -n $line | grep -v primertrim | cut -f 1 -d ":" | head -n 1 ) 2>> $err_file
//         amp_fail=$(cut -f $bedtools_column multicov.txt | awk '{{ if ( $1 < 20 ) print $0 }}' | wc -l ) 2>> $err_file
//         if [ -z "$amp_fail" ] ; then amp_fail=0 ; fi
//
//         num_of_N=$(grep -o 'N' $line*.consensus.fa | wc -l ) 2>> $err_file
//         if [ -z "$num_of_N" ] ; then num_of_N=0 ; fi
//
//         echo "$sample,$human_reads,$degenerate,${cov_and_depth[0]},${cov_and_depth[1]},$amp_fail,$num_of_N" >> covid/summary.txt
//       fi
//
//     done < <(cat !{params.sample_file} | awk '{print $1}')
//   '''
// }
//
// process file_submission {
//   publishDir "${params.outdir}", mode: 'copy'
//   tag "$sample"
//   echo true
//   cpus 1
//
//   beforeScript 'mkdir -p covid/files_for_submission logs/submission'
//
//   input:
//   file(summary) from consensus3.collect()
//
//   output:
//
//   shell:
//   '''
//   log_file=logs/submission/submission.!{workflow.sessionId}.log
//   err_file=logs/submission/submission.!{workflow.sessionId}.err
//
//   !{submission_script} !{params.outdir}
//   '''
// }
//
// process multiqc {
//   publishDir "${params.outdir}", mode: 'copy'
//   tag "$sample"
//   echo true
//   cpus 1
//
//   beforeScript 'mkdir -p covid/multiqc logs/multiqc'
//
//   input:
//   file(file) from final_summary.collect()
//
//   output:
//     path("covid/multiqc/multiqc_data")
//     file("covid/multiqc/multiqc_report.html")
//
//   shell:
//     '''
//     log_file=logs/multiqc/multiqc.!{workflow.sessionId}.log
//     err_file=logs/multiqc/multiqc.!{workflow.sessionId}.err
//
//     date | tee -a $log_file $err_file > /dev/null
//     multiqc --version >> $log_file
//
//     multiqc -f \
//         --outdir covid/multiqc \
//         !{params.outdir}/covid
//         2>> $err_file | tee -a $log_file
//     '''
// }
//
// workflow.onComplete {
//     println("Pipeline completed at: $workflow.complete")
//     println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
// }
