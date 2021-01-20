#!/usr/bin/env nextflow

println("Currently using the Cecret workflow for use with amplicon-based Illumina hybrid library prep on MiSeq\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: v.20210120")
println("")

//# nextflow run Cecret/Cecret.nf -c Cecret/config/singularity.config
//nextflow run /home/eriny/sandbox/Cecret/Cecret.nf -c /home/eriny/sandbox/Cecret/config/UPHL.config -resume -with-dag flowchart_$(date +"%H%M%S").png
// TBA plot-ampliconstats
// plot-ampliconstats results_SAMPLEID ampliconstats.txt

params.reads = workflow.launchDir + '/Sequencing_reads/Raw'
params.single_reads = workflow.launchDir + '/Sequencing_reads/Single'
if ( params.reads == params.single_reads ) {
  println("params.reads and params.single_reads can not point to the same directory!")
  println("params.reads is set to " + params.reads)
  println("params.single_reads is set to " + params.single_reads)
  exit 1
}
params.outdir = workflow.launchDir + '/cecret'

// reference files for SARS-CoV-2 (part of the github repository)
params.reference_genome = workflow.projectDir + "/config/MN908947.3.fasta"
params.gff_file = workflow.projectDir + "/config/MN908947.3.gff"
params.primer_bed = workflow.projectDir + "/config/artic_V3_nCoV-2019.bed"

params.trimmer = 'ivar'
params.cleaner = 'seqyclean'
params.aligner  = 'bwa'

// minimap2 paramaters
params.minimap2_K = '20M' // stolen from monroe

// param that coincides with the staphb/seqyclean:1.10.09 container run with singularity
params.seqyclean_contaminant_file="/Adapters_plus_PhiX_174.fasta"
params.seqyclean_minlen = 25

// for ivar
params.ivar_quality = 20
params.ivar_frequencing_threshold = 0.6
params.ivar_minimum_read_depth = 10
params.mpileup_depth = 8000

// to toggle off processes
params.prepare_reference = true
params.bcftools_variants = false
params.fastqc = true
params.ivar_variants = true
params.samtools_stats = true
params.samtools_coverage = true
params.samtools_flagstat = true
params.samtools_ampliconstats = true
params.bedtools = true
params.nextclade = true
params.pangolin = true
params.bamsnap = false // currently doesn't work. Don't turn it on until it can do non-human refrences

// for optional contamination determination
params.kraken2 = false
params.kraken2_db = ''

// for optional route of tree generation and counting snps between samples
params.relatedness = false
params.snpdists = true
params.iqtree = true
params.max_ambiguous = '0.50'
params.outgroup = 'MN908947.3'
params.mode='GTR'

params.maxcpus = Runtime.runtime.availableProcessors()
maxcpus = params.maxcpus
println("The maximum number of CPUS used in this workflow is ${maxcpus}")
if ( maxcpus < 5 ) {
  medcpus = maxcpus
} else {
  medcpus = 5
}

// this sample file contains metadata for renaming files . See README for more information
params.sample_file = workflow.launchDir + '/covid_samples.txt'
sample_file = file(params.sample_file)
params.year = Calendar.getInstance().get(Calendar.YEAR)
params.country = 'USA'
samples = []
if (sample_file.exists()) {
  println("List of COVID19 samples: " + params.sample_file)
  sample_file
    .readLines()
    .each { samples << it.split('\t')[0] }
  }

// remember to include primer files in paramaters!!!
if (file(params.primer_bed).exists()) {
  println("Primer Bed file : " + params.primer_bed)
  }
  else {
    println("A bedfile for primers is required. Set with --primer_bed or in the config file")
    exit 1
  }

// This is where the results will be
println("The files and directory for results is " + params.outdir)
println("A table summarizing results will be created: ${params.outdir}/summary.txt and ${workflow.launchDir}/run_results.txt\n")

process prepare_reference {
  publishDir "${params.outdir}", mode: 'copy', pattern: "logs/prepare_reference/*{log,err}"
  tag "reference"
  echo false
  cpus 1

  beforeScript 'mkdir -p logs/prepare_reference reference_genome'

  when:
  params.prepare_reference

  input:
  params.reference_genome

  output:
  file("reference_genome/reference.fasta*") into reference_genome
  file("logs/prepare_reference/${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/prepare_reference/!{workflow.sessionId}.log
    err_file=logs/prepare_reference/!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null

    cp !{params.reference_genome} reference_genome/reference.fasta
    bwa index reference_genome/reference.fasta
  '''
}
//reference_genome.view()

Channel
  .fromFilePairs(["${params.reads}/*_R{1,2}*.fastq.gz",
                  "${params.reads}/*_{1,2}.fastq*"], size: 2 )
  .map{ reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1], "paired" ) }
  .set { paired_reads }
//paired_reads.view()

Channel
  .fromFilePairs("${params.single_reads}/*.fastq*", size: 1 )
  .map{ reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1], "single" ) }
  .set { single_reads }

paired_reads
  .concat(single_reads)
  .ifEmpty{
    println("No fastq or fastq.gz files were found. at ${params.single_reads}")
    println("Set '--reads' to directory with paired-end reads")
    println("Set '--single_reads' to directory with single-end reads")
    exit 1
  }
  .into { fastq_reads ; fastq_reads2 ; fastq_reads3 ; fastq_reads4 }

process seqyclean {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  beforeScript 'mkdir -p seqyclean logs/seqyclean'

  when:
  params.cleaner == 'seqyclean'

  input:
  set val(sample), file(reads), val(paired_single) from fastq_reads

  output:
  tuple sample, file("seqyclean/${sample}_clean_PE{1,2}.fastq") optional true into seqyclean_paired_files
  tuple sample, file("seqyclean/${sample}_cln_SE.fastq") optional true into seqyclean_single_file
  tuple sample, file("seqyclean/${sample}_clean_PE{1,2}.fastq"), val(paired_single) optional true into seqyclean_paired_files_classification
  tuple sample, file("seqyclean/${sample}_cln_SE.fastq"), val(paired_single) optional true into seqyclean_single_file_classification
  file("seqyclean/${sample}_cl*n_SummaryStatistics.{txt,tsv}")
  file("logs/seqyclean/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(perc_kept) into seqyclean_perc_kept_results
  tuple sample, env(kept) into seqyclean_pairskept_results

  shell:
  '''
    log_file=logs/seqyclean/!{sample}.!{workflow.sessionId}.log
    err_file=logs/seqyclean/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "seqyclean version: $(seqyclean -h | grep Version)" >> $log_file

    kept=''
    perc_kept=''

    if [ "!{paired_single}" == "single" ]
    then
      seqyclean -minlen !{params.seqyclean_minlen} -qual -c !{params.seqyclean_contaminant_file} -U !{reads} -o seqyclean/!{sample}_cln 2>> $err_file >> $log_file
      kept=$(cut -f 36 seqyclean/!{sample}_cln_SummaryStatistics.tsv | grep -v "Kept" | head -n 1)
      perc_kept=$(cut -f 37 seqyclean/!{sample}_cln_SummaryStatistics.tsv | grep -v "Kept" | head -n 1)
    else
      seqyclean -minlen !{params.seqyclean_minlen} -qual -c !{params.seqyclean_contaminant_file} -1 !{reads[0]} -2 !{reads[1]} -o seqyclean/!{sample}_clean 2>> $err_file >> $log_file
      kept=$(cut -f 58 seqyclean/!{sample}_clean_SummaryStatistics.tsv | grep -v "Kept" | head -n 1)
      perc_kept=$(cut -f 59 seqyclean/!{sample}_clean_SummaryStatistics.tsv | grep -v "Kept" | head -n 1)
    fi

    if [ -z "$kept" ] ; then kept="0" ; fi
    if [ -z "$perc_kept" ] ; then perc_kept="0" ; fi
  '''
}

process fastp {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  beforeScript 'mkdir -p fastp logs/fastp'

  when:
  params.cleaner == 'fastp'

  input:
  set val(sample), file(reads), val(paired_single) from fastq_reads4

  output:
  tuple sample, file("fastp/${sample}_clean_PE{1,2}.fastq.gz") optional true into fastp_paired_files
  tuple sample, file("fastp/${sample}_cln.fastq.gz") optional true into fastp_single_file
  tuple sample, file("fastp/${sample}_clean_PE{1,2}.fastq.gz"), val(paired_single) optional true into fastp_paired_files_classification
  tuple sample, file("fastp/${sample}_cln.fastq.gz"), val(paired_single) optional true into fastp_single_file_classification
  file("fastp/${sample}_fastp.{html,json}")
  file("logs/fastp/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(passed_reads) into fastp_results

  shell:
  '''
    log_file=logs/fastp/!{sample}.!{workflow.sessionId}.log
    err_file=logs/fastp/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    fastp --version >> $log_file

    if [ "!{paired_single}" == "single" ]
    then
      fastp -i !{reads} \
        -o fastp/!{sample}_cln.fastq.gz \
        -h fastp/!{sample}_fastp.html \
        -j fastp/!{sample}_fastp.json \
        2>> $err_file >> $log_file
    else
      fastp -i !{reads[0]} \
        -I !{reads[1]} \
        -o fastp/!{sample}_clean_PE1.fastq.gz \
        -O fastp/!{sample}_clean_PE2.fastq.gz \
        -h fastp/!{sample}_fastp.html \
        -j fastp/!{sample}_fastp.json \
        2>> $err_file >> $log_file
    fi

    passed_reads=$(grep "reads passed filter" $err_file | tail -n 1 | cut -f 2 -d ":" | sed 's/ //g' )
    if [ -z "$passed_reads" ] ; then passed_reads="0" ; fi
  '''
}

seqyclean_paired_files
  .concat(fastp_paired_files)
  .concat(seqyclean_single_file)
  .concat(fastp_single_file)
  .into { clean_reads ; clean_reads2 }

seqyclean_paired_files_classification
  .concat(fastp_paired_files_classification)
  .concat(seqyclean_single_file_classification)
  .concat(fastp_single_file_classification)
  .set { clean_read_classification }

process bwa {
  publishDir "${params.outdir}", mode: 'copy', pattern: "logs/bwa/*.{log,err}"
  tag "${sample}"
  echo false
  cpus maxcpus

  beforeScript 'mkdir -p aligned logs/bwa'

  when:
  params.aligner == 'bwa'

  input:
  set val(sample), file(reads) from clean_reads
  file(reference_genome) from reference_genome

  output:
  tuple sample, file("aligned/${sample}.sam") into bwa_sams
  file("logs/bwa/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(bwa_version) into bwa_version

  shell:
  '''
    log_file=logs/bwa/!{sample}.!{workflow.sessionId}.log
    err_file=logs/bwa/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "bwa $(bwa 2>&1 | grep Version )" >> $log_file
    bwa_version="bwa : "$(bwa 2>&1 | grep Version)

    # bwa mem command
    bwa mem -t !{task.cpus} reference.fasta !{reads} 2>> $err_file > aligned/!{sample}.sam
  '''
}
//bams3.view()

process minimap2 {
  publishDir "${params.outdir}", mode: 'copy', pattern: "logs/minimap2/*.{log,err}"
  tag "${sample}"
  echo false
  cpus maxcpus

  beforeScript 'mkdir -p aligned logs/minimap2'

  when:
  params.aligner == 'minimap2'

  input:
  set val(sample), file(reads) from clean_reads2
  file(reference_genome) from reference_genome

  output:
  tuple sample, file("aligned/${sample}.sam") into minimap2_sams
  file("logs/minimap2/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(minimap2_version) into minimap2_version

  shell:
  '''
    log_file=logs/minimap2/!{sample}.!{workflow.sessionId}.log
    err_file=logs/minimap2/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    minimap2 --version >> $log_file
    minimap2_version=$(echo "minimap2 : "$(minimap2 --version))

    minimap2 -K !{params.minimap2_K} -ax sr -t !{task.cpus} -o aligned/!{sample}.sam reference.fasta !{reads} 2>> $err_file >> $log_file
  '''
}
//bams3.view()

bwa_version
  .concat(minimap2_version)
  .set { aligner_version }

bwa_sams
  .concat(minimap2_sams)
  .set { sams }

process fastqc {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo false
  cpus 1

  beforeScript 'mkdir -p fastqc logs/fastqc'

  when:
  params.fastqc

  input:
  set val(sample), file(raw), val(type) from fastq_reads2

  output:
  file("fastqc/*.{html,zip}")
  tuple sample, env(raw_1) into fastqc_1_results
  tuple sample, env(raw_2) into fastqc_2_results
  file("logs/fastqc/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/fastqc/!{sample}.!{workflow.sessionId}.log
    err_file=logs/fastqc/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    fastqc --version >> $log_file

    fastqc --outdir fastqc --threads !{task.cpus} !{raw} 2>> $err_file >> $log_file

    zipped_fastq=($(ls fastqc/*fastqc.zip) "")

    raw_1=$(unzip -p ${zipped_fastq[0]} */fastqc_data.txt | grep "Total Sequences" | awk '{ print $3 }' )
    raw_2=NA
    if [ -f "${zipped_fastq[1]}" ] ; then raw_2=$(unzip -p fastqc/*fastqc.zip */fastqc_data.txt | grep "Total Sequences" | awk '{ print $3 }' ) ; fi

    if [ -z "$raw_1" ] ; then raw_1="0" ; fi
    if [ -z "$raw_2" ] ; then raw_2="0" ; fi
  '''
}

process sort {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus maxcpus

  beforeScript 'mkdir -p aligned logs/sort'

  input:
  set val(sample), file(sam) from sams

  output:
  tuple sample, file("aligned/${sample}.sorted.bam") into pre_trim_bams, pre_trim_bams2, pre_trim_bams3
  file("logs/sort/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/sort/!{sample}.!{workflow.sessionId}.log
    err_file=logs/sort/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    samtools sort !{sam} 2>> $err_file | \
      samtools view -F 4 -o aligned/!{sample}.sorted.bam 2>> $err_file >> $log_file

    # indexing the bams
    samtools index aligned/!{sample}.sorted.bam 2>> $err_file >> $log_file
  '''
}
//bams3.view()

process ivar_trim {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  beforeScript 'mkdir -p ivar_trim logs/ivar_trim'

  when:
  params.trimmer == 'ivar'

  input:
  set val(sample), file(bam) from pre_trim_bams

  output:
  tuple sample, file("ivar_trim/${sample}.primertrim.sorted.bam") into ivar_bams
  tuple sample, file("ivar_trim/${sample}.primertrim.sorted.bam"), file("ivar_trim/${sample}.primertrim.sorted.bam.bai") into ivar_bam_bai
  file("logs/ivar_trim/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/ivar_trim/!{sample}.!{workflow.sessionId}.log
    err_file=logs/ivar_trim/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    ivar version >> $log_file

    # trimming the reads
    ivar trim -e -i !{bam} -b !{params.primer_bed} -p ivar_trim/!{sample}.primertrim 2>> $err_file >> $log_file

    # sorting and indexing the trimmed bams
    samtools sort ivar_trim/!{sample}.primertrim.bam -o ivar_trim/!{sample}.primertrim.sorted.bam 2>> $err_file >> $log_file
    samtools index ivar_trim/!{sample}.primertrim.sorted.bam 2>> $err_file >> $log_file
  '''
}

process samtools_trim {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  beforeScript 'mkdir -p samtools_trim logs/samtools_trim'

  when:
  params.trimmer == 'samtools'

  input:
  set val(sample), file(bam) from pre_trim_bams2

  output:
  tuple sample, file("samtools_trim/${sample}.primertrim.sorted.bam") into samtools_bams
  tuple sample, file("samtools_trim/${sample}.primertrim.sorted.bam"), file("samtools_trim/${sample}.primertrim.sorted.bam.bai") into samtools_bam_bai
  file("logs/samtools_trim/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/samtools_trim/!{sample}.!{workflow.sessionId}.log
    err_file=logs/samtools_trim/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    # trimming the reads
    samtools ampliconclip -b !{params.primer_bed} !{bam} 2>> $err_file | \
      samtools sort 2>> $err_file |  \
      samtools view -F 4 -o samtools_trim/!{sample}.primertrim.sorted.bam 2>> $err_file >> $log_file

    samtools index samtools_trim/!{sample}.primertrim.sorted.bam 2>> $err_file >> $log_file
  '''
}
//samtools_bams.view()

ivar_bams
  .concat(samtools_bams)
  .into { trimmed_bams ; trimmed_bams2 ; trimmed_bams3 ; trimmed_bams4 ; trimmed_bams5 }
//trimmed_bams3.view()

ivar_bam_bai
  .concat(samtools_bam_bai)
  .set { trimmed_bam_bai }
//trimmed_bam_bai2.view()

process ivar_variants {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  beforeScript 'mkdir -p ivar_variants logs/ivar_variants'

  when:
  params.ivar_variants

  input:
  set val(sample), file(bam) from trimmed_bams

  output:
  tuple sample, bam, file("ivar_variants/${sample}.variants.tsv") into ivar_variant_file
  file("logs/ivar_variants/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(variants_num) into ivar_variants_results

  shell:
  '''
    log_file=logs/ivar_variants/!{sample}.!{workflow.sessionId}.log
    err_file=logs/ivar_variants/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file
    ivar version >> $log_file

    samtools mpileup -A -d !{params.mpileup_depth} -B -Q 0 --reference !{params.reference_genome} !{bam} 2>> $err_file | \
      ivar variants -p ivar_variants/!{sample}.variants -q !{params.ivar_quality} -t !{params.ivar_frequencing_threshold} -m !{params.ivar_minimum_read_depth} -r !{params.reference_genome} -g !{params.gff_file} 2>> $err_file >> $log_file

    variants_num=$(grep "TRUE" ivar_variants/!{sample}.variants.tsv | wc -l)

    if [ -z "$variants_num" ] ; then variants_num="0" ; fi
  '''
}

process ivar_consensus {
  publishDir "${params.outdir}", mode: 'copy', pattern: "logs/ivar_consensus/*.{log,err}"
  publishDir "${params.outdir}", mode: 'copy', pattern: "consensus/*.consensus.fa"
  tag "${sample}"
  echo false
  cpus 1

  beforeScript 'mkdir -p consensus/qc_consensus/{15000,25000} logs/ivar_consensus'

  input:
  set val(sample), file(bam) from trimmed_bams2
  params.reference_genome

  output:
  tuple sample, file("consensus/${sample}.consensus.fa") into consensus, consensus2
  tuple sample, file("consensus/qc_consensus/15000/${sample}.consensus.fa") optional true into qc_consensus_15000, qc_consensus_15000_mafft
  tuple sample, file("consensus/qc_consensus/25000/${sample}.consensus.fa") optional true into qc_consensus_25000
  file("logs/ivar_consensus/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(num_N), env(num_ACTG), env(num_degenerate), env(num_total) into consensus_results
  tuple sample, env(ivar_version) into ivar_version

  shell:
  '''
    log_file=logs/ivar_consensus/!{sample}.!{workflow.sessionId}.log
    err_file=logs/ivar_consensus/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file
    ivar version >> $log_file
    ivar_version=$(ivar version | grep "version")

    samtools mpileup -A -d !{params.mpileup_depth} -B -Q 0 --reference !{params.reference_genome} !{bam} 2>> $err_file | \
      ivar consensus -q !{params.ivar_quality} -t !{params.ivar_frequencing_threshold} -m !{params.ivar_minimum_read_depth} -p consensus/!{sample}.consensus -n N 2>> $err_file >> $log_file

    num_N=$(grep -v ">" consensus/!{sample}.consensus.fa | grep -o 'N' | wc -l )
    if [ -z "$num_N" ] ; then num_N="0" ; fi

    num_ACTG=$(grep -v ">" consensus/!{sample}.consensus.fa | grep -o -E "C|A|T|G" | wc -l )
    if [ -z "$num_ACTG" ] ; then num_ACTG="0" ; fi
    if [ "$num_ACTG" -gt 15000 ] ; then cp consensus/!{sample}.consensus.fa consensus/qc_consensus/15000/!{sample}.consensus.fa ; fi
    if [ "$num_ACTG" -gt 25000 ] ; then cp consensus/!{sample}.consensus.fa consensus/qc_consensus/25000/!{sample}.consensus.fa ; fi

    num_degenerate=$(grep -v ">" consensus/!{sample}.consensus.fa | grep -o -E "B|D|E|F|H|I|J|K|L|M|O|P|Q|R|S|U|V|W|X|Y|Z" | wc -l )
    if [ -z "$num_degenerate" ] ; then num_degenerate="0" ; fi

    num_total=$(( $num_N + $num_degenerate + $num_ACTG ))
  '''
}

process bamsnap {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  beforeScript "mkdir -p bamsnap/${sample} logs/bamsnap"

  when:
  params.bamsnap

  input:
  tuple val(sample), file(bam), file(variants) from ivar_variant_file

  output:
  tuple sample, file("bamsnap/*png")
  file("logs/bamsnap/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/bamsnap/!{sample}.!{workflow.sessionId}.log
    err_file=logs/bamsnap/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    bamsnap --version >> $log_file

    bamsnap_variants=($(grep -v REGION !{variants} | awk '{ print $1 ":" $2 }' ))

    bamsnap -bam !{bam} -pos ${bamsnap_variants[@]} -out bamsnap/!{sample}/variant.png
  '''
}

process bcftools_variants {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  beforeScript 'mkdir -p bcftools_variants logs/bcftools_variants'

  when:
  params.bcftools_variants

  input:
  set val(sample), file(bam) from trimmed_bams3

  output:
  file("bcftools_variants/${sample}.vcf")
  file("logs/bcftools_variants/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(variants_num) into bcftools_variants_results

  shell:
  '''
    log_file=logs/bcftools_variants/!{sample}.!{workflow.sessionId}.log
    err_file=logs/bcftools_variants/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    bcftools --version >> $log_file

    bcftools mpileup -A -d !{params.mpileup_depth} -B -Q 0 -f !{params.reference_genome} !{bam} 2>> $err_file | \
      bcftools call -mv -Ov -o bcftools_variants/!{sample}.vcf 2>> $err_file >> $log_file

    variants_num=$(grep -v "#" bcftools_variants/!{sample}.vcf | wc -l)
    if [ -z "$variants_num" ] ; then variants_num="0" ; fi
  '''
}
//bcftools_variants.view()

pre_trim_bams3
   .combine(trimmed_bams4, by: 0)
   .into { pre_post_bams ; pre_post_bams2 ; pre_post_bams3 }
//pre_post_bams2.view()

process samtools_stats {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  beforeScript 'mkdir -p samtools_stats/aligned samtools_stats/trimmed logs/samtools_stats'

  when:
  params.samtools_stats

  input:
  set val(sample), file(aligned), file(trimmed) from pre_post_bams

  output:
  file("samtools_stats/aligned/${sample}.stats.txt")
  file("samtools_stats/trimmed/${sample}.stats.trim.txt")
  file("logs/samtools_stats/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/samtools_stats/!{sample}.!{workflow.sessionId}.log
    err_file=logs/samtools_stats/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    samtools stats !{aligned} 2>> $err_file > samtools_stats/aligned/!{sample}.stats.txt
    samtools stats !{trimmed} 2>> $err_file > samtools_stats/trimmed/!{sample}.stats.trim.txt
  '''
}

process samtools_coverage {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  beforeScript 'mkdir -p samtools_coverage/aligned samtools_coverage/trimmed logs/samtools_coverage'

  when:
  params.samtools_coverage

  input:
  set val(sample), file(aligned), file(trimmed) from pre_post_bams2

  output:
  file("samtools_coverage/aligned/${sample}.cov.{txt,hist}")
  file("samtools_coverage/trimmed/${sample}.cov.trim.{txt,hist}")
  file("logs/samtools_coverage/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(coverage) into samtools_coverage_results
  tuple sample, env(depth) into samtools_depth_results

  shell:
  '''
    log_file=logs/samtools_coverage/!{sample}.!{workflow.sessionId}.log
    err_file=logs/samtools_coverage/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    samtools coverage !{aligned} -m -o samtools_coverage/aligned/!{sample}.cov.hist 2>> $err_file >> $log_file
    samtools coverage !{aligned} -o samtools_coverage/aligned/!{sample}.cov.txt 2>> $err_file >> $log_file
    samtools coverage !{trimmed} -m -o samtools_coverage/trimmed/!{sample}.cov.trim.hist 2>> $err_file >> $log_file
    samtools coverage !{trimmed} -o samtools_coverage/trimmed/!{sample}.cov.trim.txt 2>> $err_file >> $log_file

    coverage=$(cut -f 6 samtools_coverage/trimmed/!{sample}.cov.trim.txt | tail -n 1)
    depth=$(cut -f 7 samtools_coverage/trimmed/!{sample}.cov.trim.txt | tail -n 1)
    if [ -z "$coverage" ] ; then coverage_trim="0" ; fi
    if [ -z "$depth" ] ; then depth_trim="0" ; fi
  '''
}

process samtools_flagstat {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  beforeScript 'mkdir -p samtools_flagstat/aligned samtools_flagstat/trimmed logs/samtools_flagstat'

  input:
  set val(sample), file(aligned), file(trimmed) from pre_post_bams3

  when:
  params.samtools_flagstat

  output:
  file("samtools_flagstat/{aligned,trimmed}/${sample}.flagstat.txt")
  file("logs/samtools_flagstat/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/samtools_flagstat/!{sample}.!{workflow.sessionId}.log
    err_file=logs/samtools_flagstat/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    samtools flagstat !{aligned} 2>> $err_file > samtools_flagstat/aligned/!{sample}.flagstat.txt
    samtools flagstat !{trimmed} 2>> $err_file > samtools_flagstat/trimmed/!{sample}.flagstat.txt
  '''
}

process kraken2 {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus maxcpus

  beforeScript 'mkdir -p kraken2 logs/kraken2'

  when:
  params.kraken2

  input:
  set val(sample), file(clean), val(paired_single) from clean_read_classification

  output:
  file("kraken2/${sample}_kraken2_report.txt")
  file("logs/kraken2/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(percentage_cov) into kraken2_sars_results
  tuple sample, env(percentage_human) into kraken2_human_results

  shell:
  '''
    log_file=logs/kraken2/!{sample}.!{workflow.sessionId}.log
    err_file=logs/kraken2/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    kraken2 --version >> $log_file

    if [ "!{paired_single}" == "single" ]
    then
      kraken2 \
        --classified-out cseqs#.fq \
        --threads !{task.cpus} \
        --db !{params.kraken2_db} \
        !{clean} \
        --report kraken2/!{sample}_kraken2_report.txt \
        2>> $err_file >> $log_file
    else
      kraken2 --paired \
        --classified-out cseqs#.fq \
        --threads !{task.cpus} \
        --db !{params.kraken2_db} \
        !{clean} \
        --report kraken2/!{sample}_kraken2_report.txt \
        2>> $err_file >> $log_file
    fi

    percentage_human=$(grep "Homo sapiens" kraken2/!{sample}_kraken2_report.txt | awk '{print $1}')
    percentage_cov=$(grep "Severe acute respiratory syndrome coronavirus 2" kraken2/!{sample}_kraken2_report.txt | awk '{print $1}')

    if [ -z "$percentage_human" ] ; then percentage_human="0" ; fi
    if [ -z "$percentage_cov" ] ; then percentage_cov="0" ; fi
  '''
}
//kraken2_results2.view()

process bedtools {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  beforeScript 'mkdir -p bedtools logs/bedtools'

  when:
  params.bedtools

  input:
  set val(sample), file(bam), file(bai) from trimmed_bam_bai
  params.primer_bed

  output:
  file("bedtools/${sample}.multicov.txt")
  tuple sample, env(num_failed_amplicons) into bedtools_results
  file("logs/bedtools/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/bedtools/!{sample}.!{workflow.sessionId}.log
    err_file=logs/bedtools/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    bedtools --version >> $log_file

    cat !{params.primer_bed} | \
      grep -v "alt" | \
      awk '{ if ($0 ~ "LEFT") { print $1 "\t" $2 } else {print $3 "\t" $4 "\t" $5 }}' | \
      paste - - | \
      sed 's/_RIGHT//g' > amplicon.bed

    bedtools multicov -bams !{bam} -bed amplicon.bed 2>> $err_file >> bedtools/!{sample}.multicov.txt

    num_failed_amplicons=$(cut -f 6 bedtools/!{sample}.multicov.txt | awk '{ if ( $1 < 20 ) print $0 }' | wc -l )
    if [ -z "$num_failed_amplicons" ] ; then num_failed_amplicons=0 ; fi
  '''
}

process samtools_ampliconstats {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  beforeScript 'mkdir -p samtools_ampliconstats logs/samtools_ampliconstats'

  when:
  params.samtools_ampliconstats

  input:
  set val(sample), file(bam) from trimmed_bams5
  params.primer_bed

  output:
  file("samtools_ampliconstats/${sample}_ampliconstats.txt")
  tuple sample, env(num_failed_amplicons) into samtools_ampliconstats_results
  file("logs/samtools_ampliconstats/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/samtools_ampliconstats/!{sample}.!{workflow.sessionId}.log
    err_file=logs/samtools_ampliconstats/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    samtools ampliconstats !{params.primer_bed} !{bam} 2>> $err_file > samtools_ampliconstats/!{sample}_ampliconstats.txt

    num_failed_amplicons=$(grep ^FREADS samtools_ampliconstats/!{sample}_ampliconstats.txt | cut -f 2- | tr '\t' '\n' | awk '{ if ($1 < 20) print $0 }' | wc -l)
    if [ -z "$num_failed_amplicons" ] ; then num_failed_amplicons=0 ; fi
  '''
}

process pangolin {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus medcpus

  beforeScript 'mkdir -p pangolin logs/pangolin'

  when:
  params.pangolin

  input:
  set val(sample), file(fasta) from consensus

  output:
  file("pangolin/${sample}/lineage_report.csv")
  tuple sample, env(pangolin_lineage) into pangolin_lineage_results
  tuple sample, env(pangolin_status) into pangolin_status_results
  file("logs/pangolin/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/pangolin/!{sample}.!{workflow.sessionId}.log
    err_file=logs/pangolin/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    pangolin --version >> $log_file
    pangolin -lv >> $log_file

    pangolin --threads !{task.cpus} --outdir pangolin/!{sample} !{fasta} 2>> $err_file >> $log_file

    pangolin_lineage=$(tail -n 1 pangolin/!{sample}/lineage_report.csv | cut -f 2 -d "," | grep -v "lineage" )
    pangolin_status=$(tail -n 1 pangolin/!{sample}/lineage_report.csv | cut -f 5 -d "," )

    if [ -z "$pangolin_lineage" ] ; then pangolin_lineage="NA" ; fi
    if [ -z "$pangolin_status" ] ; then pangolin_status="NA" ; fi
  '''
}

process nextclade {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus medcpus

  beforeScript 'mkdir -p nextclade logs/nextclade'

  when:
  params.nextclade

  input:
  set val(sample), file(fasta) from consensus2

  output:
  file("nextclade/${sample}_nextclade_report.csv")
  tuple sample, env(nextclade_clade) into nextclade_clade_results
  file("logs/nextclade/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/nextclade/!{sample}.!{workflow.sessionId}.log
    err_file=logs/nextclade/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    nextclade --version >> $log_file

    nextclade --jobs !{task.cpus} --input-fasta !{fasta} --output-csv nextclade/!{sample}_nextclade_report.csv 2>> $err_file >> $log_file
    nextclade_clade=$(cat nextclade/!{sample}_nextclade_report.csv | grep !{sample} | cut -f 2 -d ";" | head -n 1 )
    if [ -z "$nextclade_clade" ] ; then nextclade_clade="clade" ; fi
  '''
}

consensus_results
//tuple sample, env(num_N), env(num_ACTG), env(num_degenerate), env(num_total) into consensus_results
  .join(fastqc_1_results, remainder: true, by: 0)
  .join(fastqc_2_results, remainder: true, by: 0)
  .join(seqyclean_pairskept_results, remainder: true, by: 0)
  .join(seqyclean_perc_kept_results, remainder: true, by: 0)
  .join(fastp_results, remainder: true, by: 0)
  .join(ivar_variants_results, remainder: true, by: 0)
  .join(bcftools_variants_results, remainder: true, by:0)
  .join(samtools_coverage_results, remainder: true, by: 0)
  .join(samtools_depth_results, remainder: true, by: 0)
  .join(kraken2_human_results, remainder: true, by: 0)
  .join(kraken2_sars_results, remainder: true, by: 0)
  .join(pangolin_lineage_results, remainder: true, by: 0)
  .join(pangolin_status_results, remainder: true, by: 0)
  .join(nextclade_clade_results, remainder: true, by: 0)
  .join(bedtools_results, remainder: true, by: 0)
  .join(samtools_ampliconstats_results, remainder: true, by: 0)
  .join(aligner_version, remainder: true, by:0)
  .join(ivar_version, remainder: true, by: 0)
  .set { results }
//results2.view()

process summary {
  publishDir "${params.outdir}", mode: 'copy', overwrite: true
  tag "${sample}"
  echo false
  cpus 1

  beforeScript 'mkdir -p summary logs/summary'

  input:
  set val(sample), val(num_N), val(num_ACTG), val(num_degenerate), val(num_total),
    val(raw_1),
    val(raw_2),
    val(pairskept),
    val(perc_kept),
    val(reads_passed),
    val(ivar_variants),
    val(bcftools_variants),
    val(coverage),
    val(depth),
    val(percentage_human),
    val(percentage_cov),
    val(pangolin_lineage),
    val(pangolin_status),
    val(nextclade_clade),
    val(bedtools_num_failed_amplicons),
    val(samtools_num_failed_amplicons),
    val(bwa_version),
    val(ivar_version) from results

  output:
  file("summary/${sample}.summary.txt") into summary
  file("logs/summary/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/summary/!{sample}.!{workflow.sessionId}.log
    err_file=logs/summary/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null

    sample_id=$(echo !{sample} | cut -f 1 -d "-" )

    echo -e "sample_id\tsample\taligner_version\tivar_version\tpangolin_lineage\tpangolin_status\tnextclade_clade\tfastqc_raw_reads_1\tfastqc_raw_reads_2\tseqyclean_pairs_kept_after_cleaning\tseqyclean_percent_kept_after_cleaning\tfastp_reads_passed\tdepth_after_trimming\tcoverage_after_trimming\t%_human_reads\t%_SARS-COV-2_reads\tivar_num_variants_identified\tbcftools_variants_identified\tbedtools_num_failed_amplicons\tsamtools_num_failed_amplicons\tnum_N\tnum_degenerage\tnum_ACTG\tnum_total" > summary/!{sample}.summary.txt
    echo -e "${sample_id}\t!{sample}\t!{bwa_version}\t!{ivar_version}\t!{pangolin_lineage}\t!{pangolin_status}\t!{nextclade_clade}\t!{raw_1}\t!{raw_2}\t!{pairskept}\t!{perc_kept}\t!{reads_passed}\t!{depth}\t!{coverage}\t!{percentage_human}\t!{percentage_cov}\t!{ivar_variants}\t!{bcftools_variants}\t!{bedtools_num_failed_amplicons}\t!{samtools_num_failed_amplicons}\t!{num_N}\t!{num_degenerate}\t!{num_ACTG}\t!{num_total}" >> summary/!{sample}.summary.txt
  '''
}

process combined_summary {
  publishDir "${params.outdir}", mode: 'copy', overwrite: true, pattern: "summary.txt"
  publishDir "${params.outdir}", mode: 'copy', overwrite: true, pattern: "logs/summary/*.{log,err}"
  publishDir "${workflow.launchDir}", mode: 'copy', overwrite: true, pattern: "run_results.txt"
  tag "summary"
  echo false
  cpus 1

  beforeScript 'mkdir -p submission_files logs/summary'

  input:
  file(summary) from summary.collect()

  output:
  file("summary.txt")
  file("run_results.txt")
  file("logs/summary/summary.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/summary/summary.!{workflow.sessionId}.log
    err_file=logs/summary/summary.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null

    cat *.summary.txt | grep "sample_id" | head -n 1 > summary.txt
    cat *summary.txt | grep -v "sample_id" | sort | uniq >> summary.txt 2>> $err_file

    cp summary.txt run_results.txt
  '''
}

process mafft {
  publishDir "${params.outdir}", mode: 'copy'
  tag "Multiple Sequence Alignment"
  echo false
  cpus maxcpus

  beforeScript 'mkdir -p mafft logs/mafft'

  input:
  file(consensus) from qc_consensus_15000_mafft.collect()
  params.reference_genome

  output:
  file("mafft/mafft_aligned.fasta") into msa_file
  file("mafft/mafft_aligned.fasta") into msa_file2
  file("logs/mafft/mafft.${workflow.sessionId}.{log,err}")

  when:
  params.relatedness

  shell:
  '''
    log_file=logs/mafft/mafft.!{workflow.sessionId}.log
    err_file=logs/mafft/mafft.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    echo "mafft version:" >> $log_file
    mafft --version 2>&1 >> $log_file

    echo ">!{params.outgroup}" > reference.fasta
    grep -v ">" !{params.reference_genome} >> reference.fasta

    cat *fa > ultimate_consensus.fasta
    mafft --auto \
      --thread !{task.cpus} \
      --maxambiguous !{params.max_ambiguous} \
      --addfragments ultimate_consensus.fasta \
      reference.fasta \
      > mafft/mafft_aligned.fasta \
      2>> $err_file
  '''
}

process snpdists {
  publishDir "${params.outdir}", mode: 'copy'
  tag "snp-dists"
  echo false
  cpus medcpus

  beforeScript 'mkdir -p snp-dists logs/snp-dists'

  input:
  file(msa) from msa_file

  output:
  file("snp-dists/snp-dists.txt")
  file("logs/snp-dists/snp-dists.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/snp-dists/snp-dists.!{workflow.sessionId}.log
    err_file=logs/snp-dists/snp-dists.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    snp-dists -v >> $log_file

    snp-dists !{msa} > snp-dists/snp-dists.txt 2> $err_file
  '''
}
//msa_file2.view()

process iqtree {
  publishDir "${params.outdir}", mode: 'copy'
  tag "iqtree"
  echo false
  cpus maxcpus

  beforeScript 'mkdir -p iqtree logs/iqtree'

  input:
  file(msa) from msa_file2

  output:
  file("iqtree/iqtree.{iqtree,treefile,mldist,log}")
  file("logs/iqtree/iqtree.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/iqtree/iqtree.!{workflow.sessionId}.log
    err_file=logs/iqtree/iqtree.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    iqtree --version >> $log_file

    cat !{msa} | sed 's/!{params.outgroup}.*/!{params.outgroup}/g' > !{msa}.tmp
    mv !{msa}.tmp !{msa}

    # creating a tree
  	iqtree -ninit 2 \
      -n 2 \
      -me 0.05 \
      -nt AUTO \
      -ntmax !{task.cpus} \
      -s !{msa} \
      -pre iqtree/iqtree \
      -m !{params.mode} \
      -o !{params.outgroup} \
      >> $log_file 2>> $err_file
  '''
}

process rename_fastq {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  beforeScript 'mkdir -p submission_files logs/rename_fastq'

  input:
  set val(sample), file(reads), val(paired_single) from fastq_reads3
  sample_file

  when:
  for(int i =0; i < samples.size(); i++) {
    if(sample.contains(samples[i])) { return true }
  }

  output:
  tuple val(sample), env(sample_id), env(submission_id) into submission_ids, submission_ids2
  file("submission_files/*.fastq.gz")
  file("logs/rename_fastq/${sample}.${workflow.sessionId}.{err,log}")

  shell:
  '''
    log_file=logs/rename_fastq/!{sample}.!{workflow.sessionId}.log
    err_file=logs/rename_fastq/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null

    sample_id=''
    submission_id=''
    while read line
    do
      lab_accession=$(echo $line | awk '{print $1}' )
      if [[ "!{sample}" == *"$lab_accession"* ]]
      then
        sample_id=$lab_accession
        submission_id=$(echo $line | awk '{print $2}' )
      fi
    done < !{sample_file}

    echo !{sample}
    echo $submission_id

    if [ -z "$sample_id" ] ; then sample_id=!{sample} ; fi
    if [ -z "$submission_id" ] ; then submission_id=!{sample} ; fi

    if [ "!{paired_single}" == "single" ]
    then
      echo "single"
      cp !{reads[0]} submission_files/$submission_id.fastq.gz
    else
      cp !{reads[0]} submission_files/$submission_id.R1.fastq.gz  2>> $err_file
      cp !{reads[1]} submission_files/$submission_id.R2.fastq.gz  2>> $err_file
    fi
  '''
}
//submission_ids.view()

qc_consensus_15000
  .join(submission_ids)
  .set { ids_genbank }
//ids_genbank.view()

qc_consensus_25000
  .join(submission_ids2, by:0)
  .set { ids_gisiad }
//ids_gisiad.view()

process prepare_gisaid {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  beforeScript 'mkdir -p submission_files logs/prepare_gisaid'

  input:
  set val(sample), file(consensus), val(sample_id), val(submission_id) from ids_gisiad

  when:
  for(int i =0; i < samples.size(); i++) {
    if(sample.contains(samples[i])) { return true }
  }

  output:
  file("submission_files/${submission_id}.gisaid.fa") into gisaid_fasta
  file("logs/prepare_gisaid/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/prepare_gisaid/!{sample}.!{workflow.sessionId}.log
    err_file=logs/prepare_gisaid/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null

    # getting the consensus fasta file
    # changing the fasta header

    echo ">hCoV-19/!{params.country}/!{submission_id}/!{params.year}" > submission_files/!{submission_id}.gisaid.fa
    grep -v ">" !{consensus} | fold -w 75 >> submission_files/!{submission_id}.gisaid.fa  2>> $err_file
  '''
}

process prepare_genbank {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo false
  cpus 1

  beforeScript 'mkdir -p submission_files logs/prepare_genbank'

  input:
  set val(sample), file(consensus), val(sample_id), val(submission_id) from ids_genbank
  sample_file

  when:
  for(int i =0; i < samples.size(); i++) {
    if(sample.contains(samples[i])) { return true }
  }

  output:
  file("submission_files/${submission_id}.genbank.fa") into genbank_fasta
  file("logs/prepare_genbank/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/prepare_genbank/!{sample}.!{workflow.sessionId}.log
    err_file=logs/prepare_genbank/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null

    !{workflow.projectDir}/bin/genbank_submission.sh -f !{consensus} -m !{sample_file} -y !{params.year} -o submission_files/!{submission_id}.genbank.fa 2>> $err_file >> $log_file
  '''
}

process combine_gisaid {
  publishDir "${params.outdir}", mode: 'copy'
  tag "multifasta"
  echo false
  cpus 1

  beforeScript 'mkdir -p submission_files logs/gisaid_fasta'

  input:
  file(fastas) from gisaid_fasta.collect()

  when:
  if (file(params.sample_file).exists()) { return true }

  output:
  file("submission_files/*.gisaid_submission.fasta") optional true
  file("logs/gisaid_fasta/gisaid_fasta.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/gisaid_fasta/gisaid_fasta.!{workflow.sessionId}.log
    err_file=logs/gisaid_fasta/gisaid_fasta.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null

    run_id=$(echo "!{params.outdir}" | rev | cut -f 1 -d '/' | rev )
    run_id=${run_id: -6}
    if [ -z "$run_id" ] ; then run_id="Submission" ; fi
    cat *gisaid.fa > submission_files/$run_id.gisaid_submission.fasta 2>> $err_file
  '''
}

process combine_genbank {
  publishDir "${params.outdir}", mode: 'copy'
  tag "multifasta"
  echo false
  cpus 1

  beforeScript 'mkdir -p submission_files logs/combine_genbank'

  input:
  file(fastas) from genbank_fasta.collect()

  when:
  if (file(params.sample_file).exists()) { return true }

  output:
  file("submission_files/*.genbank_submission.fasta") optional true
  file("logs/combine_genbank/combine_genbank.${workflow.sessionId}.{log,err}")

  shell:
  '''
    log_file=logs/combine_genbank/combine_genbank.!{workflow.sessionId}.log
    err_file=logs/combine_genbank/combine_genbank.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null

    run_id=$(echo "!{params.outdir}" | rev | cut -f 1 -d '/' | rev )
    run_id=${run_id: -6}
    if [ -z "$run_id" ] ; then run_id="Submission" ; fi
    cat *genbank.fa > submission_files/$run_id.genbank_submission.fasta 2>> $err_file
  '''
}

workflow.onComplete {
    println("Pipeline completed at: $workflow.complete")
    println("A summary of results can be found in a tab-delimited file: ${workflow.launchDir}/run_results.txt")
    if (file(params.sample_file).exists()) { println("SRA, GenBank, and GISAID submission-ready files are located at ${workflow.launchDir}/submission_files") }
    println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
}
