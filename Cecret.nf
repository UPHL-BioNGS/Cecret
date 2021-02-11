#!/usr/bin/env nextflow

println("Currently using the Cecret workflow for use with amplicon-based Illumina hybrid library prep on MiSeq\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: v.20210211")
println("")

// TBA plot-ampliconstats
// plot-ampliconstats results_SAMPLEID ampliconstats.txt

params.reads = workflow.launchDir + '/Sequencing_reads/Raw'
params.single_reads = workflow.launchDir + '/Sequencing_reads/Single'
if ( params.reads == params.single_reads ) {
  println("'params.reads' and 'params.single_reads' cannot point to the same directory!")
  println("'params.reads' is set to " + params.reads)
  println("'params.single_reads' is set to " + params.single_reads)
  exit 1
}
params.outdir = workflow.launchDir + '/cecret'

// reference files for SARS-CoV-2 (part of the github repository)
params.reference_genome = workflow.projectDir + "/configs/MN908947.3.fasta"
params.gff_file = workflow.projectDir + "/configs/MN908947.3.gff"
params.primer_bed = workflow.projectDir + "/configs/artic_V3_nCoV-2019.bed"

params.trimmer = 'ivar'
params.cleaner = 'seqyclean'
params.aligner = 'bwa'

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
params.bcftools_variants = false // fails to download a lot
params.fastqc = true
params.ivar_variants = true
params.samtools_stats = true
params.samtools_coverage = true
params.samtools_flagstat = true
params.samtools_ampliconstats = true
params.bedtools = true
params.nextclade = true
params.pangolin = true
params.bamsnap = false // can be really slow
params.rename = false

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

// for optional renaming of files for GISAID and GenBank submissions
params.sample_file = workflow.launchDir + '/covid_samples.csv'
params.gisaid_threshold = '25000'
params.genbank_threshold = '15000'

params.maxcpus = Runtime.runtime.availableProcessors()
println("The maximum number of CPUS used in this workflow is ${params.maxcpus}")
if ( params.maxcpus < 5 ) {
  params.medcpus = params.maxcpus
} else {
  params.medcpus = 5
}

// This is where the results will be
println("The files and directory for results is " + params.outdir)
println("A table summarizing results will be created: ${params.outdir}/summary.txt and ${workflow.launchDir}/run_results.txt\n")

Channel
  .fromPath(params.reference_genome, type:'file')
  .ifEmpty{
    println("No reference genome was selected. Set with 'params.reference_genome'")
    exit 1
  }
  .view { "Reference Genome : $it"}
  .into { reference_genome ; reference_genome2 ; reference_genome_mafft ; reference_genome_bamsnap }

Channel
  .fromPath(params.gff_file, type:'file')
  .view { "GFF file for Reference Genome : $it"}
  .set { gff_file }

Channel
  .fromPath(params.primer_bed, type:'file')
  .ifEmpty{
    println("A bedfile for primers is required. Set with 'params.primer_bed'.")
    exit 1
  }
  .view { "Primer BedFile : $it"}
  .into { primer_bed ; primer_bed_bedtools ; primer_bed_ampliconstats }

Channel
  .fromFilePairs(["${params.reads}/*_R{1,2}*.fastq.gz",
                  "${params.reads}/*_{1,2}.fastq*"], size: 2 )
  .map{ reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1], "paired" ) }
  .set { paired_reads }

Channel
  .fromFilePairs("${params.single_reads}/*.fastq*", size: 1 )
  .map{ reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1], "single" ) }
  .set { single_reads }


sample_file = params.rename ? Channel.fromPath(params.sample_file, type:'file').view { "Sample File : $it" } : Channel.empty()

paired_reads
  .concat(single_reads)
  .ifEmpty{
    println("No fastq or fastq.gz files were found at ${params.reads} or ${params.single_reads}")
    println("Set 'params.reads' to directory with paired-end reads")
    println("Set 'params.single_reads' to directory with single-end reads")
    exit 1
  }
  .into { fastq_reads_seqyclean ; fastq_reads_fastp ; fastq_reads_fastqc ; fastq_reads_rename }

println("") // just for aesthetics

process seqyclean {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  when:
  params.cleaner == 'seqyclean'

  input:
  set val(sample), file(reads), val(paired_single) from fastq_reads_seqyclean

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
    mkdir -p seqyclean logs/seqyclean
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

  when:
  params.cleaner == 'fastp'

  input:
  set val(sample), file(reads), val(paired_single) from fastq_reads_fastp

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
    mkdir -p fastp logs/fastp
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
  .combine(reference_genome)
  .into { clean_reads_bwa ; clean_reads_minimap2 }

seqyclean_paired_files_classification
  .concat(fastp_paired_files_classification)
  .concat(seqyclean_single_file_classification)
  .concat(fastp_single_file_classification)
  .set { clean_reads_classification }

process bwa {
  publishDir "${params.outdir}", mode: 'copy', pattern: "logs/bwa/*.{log,err}"
  tag "${sample}"
  echo false
  cpus params.maxcpus

  when:
  params.aligner == 'bwa'

  input:
  set val(sample), file(reads), file(reference_genome) from clean_reads_bwa

  output:
  tuple sample, file("aligned/${sample}.sam") into bwa_sams
  file("logs/bwa/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(bwa_version) into bwa_version

  shell:
  '''
    mkdir -p aligned logs/bwa
    log_file=logs/bwa/!{sample}.!{workflow.sessionId}.log
    err_file=logs/bwa/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "bwa $(bwa 2>&1 | grep Version )" >> $log_file
    bwa_version="bwa : "$(bwa 2>&1 | grep Version)

    # index the reference fasta file
    bwa index !{reference_genome}

    # bwa mem command
    bwa mem -t !{task.cpus} !{reference_genome} !{reads} 2>> $err_file > aligned/!{sample}.sam
  '''
}

process minimap2 {
  publishDir "${params.outdir}", mode: 'copy', pattern: "logs/minimap2/*.{log,err}"
  tag "${sample}"
  echo false
  cpus params.maxcpus

  when:
  params.aligner == 'minimap2'

  input:
  set val(sample), file(reads), file(reference_genome) from clean_reads_minimap2

  output:
  tuple sample, file("aligned/${sample}.sam") into minimap2_sams
  file("logs/minimap2/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(minimap2_version) into minimap2_version

  shell:
  '''
    mkdir -p aligned logs/minimap2
    log_file=logs/minimap2/!{sample}.!{workflow.sessionId}.log
    err_file=logs/minimap2/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    minimap2 --version >> $log_file
    minimap2_version=$(echo "minimap2 : "$(minimap2 --version))

    minimap2 -K !{params.minimap2_K} -ax sr -t !{task.cpus} -o aligned/!{sample}.sam !{reference_genome} !{reads} 2>> $err_file >> $log_file
  '''
}

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

  when:
  params.fastqc

  input:
  set val(sample), file(raw), val(type) from fastq_reads_fastqc

  output:
  file("fastqc/*.{html,zip}")
  tuple sample, env(raw_1) into fastqc_1_results
  tuple sample, env(raw_2) into fastqc_2_results
  file("logs/fastqc/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p fastqc logs/fastqc
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
  cpus params.maxcpus

  input:
  set val(sample), file(sam) from sams

  output:
  tuple sample, file("aligned/${sample}.sorted.bam") into pre_trim_bams, pre_trim_bams2
  tuple sample, file("aligned/${sample}.sorted.bam"), file("aligned/${sample}.sorted.bam.bai") into pre_trim_bams_bamsnap
  file("logs/sort/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p aligned logs/sort
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

pre_trim_bams
  .combine(primer_bed)
  .into {pre_trim_bams_ivar ; pre_trim_bams_samtools }

process ivar_trim {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  when:
  params.trimmer == 'ivar'

  input:
  set val(sample), file(bam), file(primer_bed) from pre_trim_bams_ivar

  output:
  tuple sample, file("ivar_trim/${sample}.primertrim.sorted.bam") into ivar_bams
  tuple sample, file("ivar_trim/${sample}.primertrim.sorted.bam"), file("ivar_trim/${sample}.primertrim.sorted.bam.bai") into ivar_bam_bai
  file("logs/ivar_trim/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p ivar_trim logs/ivar_trim
    log_file=logs/ivar_trim/!{sample}.!{workflow.sessionId}.log
    err_file=logs/ivar_trim/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    ivar version >> $log_file

    # trimming the reads
    ivar trim -e -i !{bam} -b !{primer_bed} -p ivar_trim/!{sample}.primertrim 2>> $err_file >> $log_file

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

  when:
  params.trimmer == 'samtools'

  input:
  set val(sample), file(bam), file(primer_bed) from pre_trim_bams_samtools

  output:
  tuple sample, file("samtools_trim/${sample}.primertrim.sorted.bam") into samtools_bams
  tuple sample, file("samtools_trim/${sample}.primertrim.sorted.bam"), file("samtools_trim/${sample}.primertrim.sorted.bam.bai") into samtools_bam_bai
  file("logs/samtools_trim/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p samtools_trim logs/samtools_trim
    log_file=logs/samtools_trim/!{sample}.!{workflow.sessionId}.log
    err_file=logs/samtools_trim/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    # trimming the reads
    samtools ampliconclip -b !{primer_bed} !{bam} 2>> $err_file | \
      samtools sort 2>> $err_file |  \
      samtools view -F 4 -o samtools_trim/!{sample}.primertrim.sorted.bam 2>> $err_file >> $log_file

    samtools index samtools_trim/!{sample}.primertrim.sorted.bam 2>> $err_file >> $log_file
  '''
}

ivar_bams
  .concat(samtools_bams)
  .into { trimmed_bams ; trimmed_bams4 ; trimmed_bams5 }

trimmed_bams5
  .combine(primer_bed_ampliconstats)
  .set { trimmed_bams_ampliconstats }

trimmed_bams
 .combine(reference_genome2)
 .into { trimmed_bams_genome ; trimmed_bams_ivar_consensus ; trimmed_bams_bcftools_variants }

trimmed_bams_genome
 .combine(gff_file)
 .set { trimmed_bams_ivar_variants }

ivar_bam_bai
  .concat(samtools_bam_bai)
  .combine(primer_bed_bedtools)
  .set { trimmed_bam_bai }

process ivar_variants {
  publishDir "${params.outdir}", mode: 'copy',  pattern: "logs/ivar_variants/*.{log,err}"
  publishDir "${params.outdir}", mode: 'copy',  pattern: "ivar_variants/*.tsv"
  tag "${sample}"
  echo false
  cpus 1

  when:
  params.ivar_variants

  input:
  set val(sample), file(bam), file(reference_genome), file(gff_file) from trimmed_bams_ivar_variants

  output:
  tuple sample, file("ivar_variants/${sample}.variants.tsv") into ivar_variant_file
  file("logs/ivar_variants/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(variants_num) into ivar_variants_results

  shell:
  '''
    mkdir -p ivar_variants logs/ivar_variants
    log_file=logs/ivar_variants/!{sample}.!{workflow.sessionId}.log
    err_file=logs/ivar_variants/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file
    ivar version >> $log_file

    samtools mpileup -A -d !{params.mpileup_depth} -B -Q 0 --reference !{reference_genome} !{bam} 2>> $err_file | \
      ivar variants -p ivar_variants/!{sample}.variants -q !{params.ivar_quality} -t !{params.ivar_frequencing_threshold} -m !{params.ivar_minimum_read_depth} -r !{reference_genome} -g !{gff_file} 2>> $err_file >> $log_file

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

  input:
  set val(sample), file(bam), file(reference_genome) from trimmed_bams_ivar_consensus

  output:
  tuple sample, file("consensus/${sample}.consensus.fa") into consensus_pangolin, consensus_nextclade
  tuple sample, file("consensus/${sample}.consensus.fa"), env(num_ACTG) into consensus_rename
  tuple sample, file("consensus/qc_consensus/15000/${sample}.consensus.fa") optional true into qc_consensus_15000_mafft
  file("logs/ivar_consensus/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(num_N), env(num_ACTG), env(num_degenerate), env(num_total) into consensus_results
  tuple sample, env(ivar_version) into ivar_version

  shell:
  '''
    mkdir -p consensus/qc_consensus/{15000,25000} logs/ivar_consensus
    log_file=logs/ivar_consensus/!{sample}.!{workflow.sessionId}.log
    err_file=logs/ivar_consensus/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file
    ivar version >> $log_file
    ivar_version=$(ivar version | grep "version")

    samtools mpileup -A -d !{params.mpileup_depth} -B -Q 0 --reference !{reference_genome} !{bam} 2>> $err_file | \
      ivar consensus -q !{params.ivar_quality} -t !{params.ivar_frequencing_threshold} -m !{params.ivar_minimum_read_depth} -p consensus/!{sample}.consensus -n N 2>> $err_file >> $log_file

    num_N=$(grep -v ">" consensus/!{sample}.consensus.fa | grep -o 'N' | wc -l )
    if [ -z "$num_N" ] ; then num_N="0" ; fi

    num_ACTG=$(grep -v ">" consensus/!{sample}.consensus.fa | grep -o -E "C|A|T|G" | wc -l )
    if [ -z "$num_ACTG" ] ; then num_ACTG="0" ; fi
    if [ "$num_ACTG" -gt 15000 ] ; then cp consensus/!{sample}.consensus.fa consensus/qc_consensus/15000/!{sample}.consensus.fa ; fi

    num_degenerate=$(grep -v ">" consensus/!{sample}.consensus.fa | grep -o -E "B|D|E|F|H|I|J|K|L|M|O|P|Q|R|S|U|V|W|X|Y|Z" | wc -l )
    if [ -z "$num_degenerate" ] ; then num_degenerate="0" ; fi

    num_total=$(( $num_N + $num_degenerate + $num_ACTG ))
  '''
}

process bcftools_variants {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  when:
  params.bcftools_variants

  input:
  set val(sample), file(bam), file(reference_genome) from trimmed_bams_bcftools_variants

  output:
  tuple sample, file("bcftools_variants/${sample}.vcf") into bcftools_variants_file
  file("logs/bcftools_variants/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(variants_num) into bcftools_variants_results

  shell:
  '''
    mkdir -p bcftools_variants logs/bcftools_variants
    log_file=logs/bcftools_variants/!{sample}.!{workflow.sessionId}.log
    err_file=logs/bcftools_variants/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    bcftools --version >> $log_file

    bcftools mpileup -A -d !{params.mpileup_depth} -B -Q 0 -f !{reference_genome} !{bam} 2>> $err_file | \
      bcftools call -mv -Ov -o bcftools_variants/!{sample}.vcf 2>> $err_file >> $log_file

    variants_num=$(grep -v "#" bcftools_variants/!{sample}.vcf | wc -l)
    if [ -z "$variants_num" ] ; then variants_num="0" ; fi
  '''
}

pre_trim_bams_bamsnap
  .join(ivar_variant_file, remainder: true, by:0)
  .join(bcftools_variants_file, remainder: true, by:0)
  .combine(reference_genome_bamsnap)
  .set { bamsnap_files }

process bamsnap {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus params.medcpus
  errorStrategy 'ignore'

  when:
  params.bamsnap

  input:
  tuple val(sample), file(bam), file(bai), file(variant_file), file(vcf), file(reference_genome) from bamsnap_files

  output:
  file("bamsnap/${sample}/{ivar,bcftools}/*.{png,log}") optional true
  file("bamsnap/${sample}/*.{png,log}") optional true
  file("logs/bamsnap/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p logs/bamsnap
    log_file=logs/bamsnap/!{sample}.!{workflow.sessionId}.log
    err_file=logs/bamsnap/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    bamsnap --version >> $log_file


    if [[ "!{variant_file}" != *"input"* ]]
    then
      reference_length=$(grep -v ">" /home/IDGenomics_NAS/testing_cecret/20210208/work/72/7d705566c975564054eadf55f3ca95/MN908947.3.fasta | wc -m)
      if [ $reference_length -lt 1050 ] ; then echo "The reference isn't long enough for bamsnap. Set 'bamsnap.params = false' "; exit 0 ; fi
      max_position=$(( $reference_length - 1050 ))
      echo "Reminder: the min variant position is 501 and the maximum position is $max_position" | tee -a $log_file

      mkdir -p bamsnap/!{sample}
      bamsnap_variants=($(grep TRUE !{variant_file} | awk '{ if ( $2 < 501 ) $2=501 ; print $1 " " $2 }' | awk -v max=$max_position '{ if ( $2 > max ) $2=max ; print $1 ":" $2 }' | sort | uniq ))

      for refvariant in ${bamsnap_variants[@]}
      do
        variant=$(echo $refvariant | cut -f 2 -d ":" )
        bamsnap -draw coordinates bamplot coverage base \
          -process !{task.cpus} \
          -ref !{reference_genome} \
          -bam !{bam} \
          -out bamsnap/!{sample}/ivar/$variant.png \
          -pos $refvariant \
          -imagetype png \
          -save_image_only 2>> $err_file | tee -a $log_file
      done
    fi

    if [[ "!{vcf}" != *"input"* ]]
    then
      mkdir -p bamsnap/!{sample}
      bamsnap -draw coordinates bamplot coverage base \
        -process !{task.cpus} \
        -ref !{reference_genome} \
        -bam !{bam} \
        -vcf !{vcf} \
        -out bamsnap/!{sample}/bcftools \
        -imagetype png \
        -save_image_only 2>> $err_file >> $log_file
    fi
  '''
}

pre_trim_bams2
   .combine(trimmed_bams4, by: 0)
   .into { pre_post_bams ; pre_post_bams2 ; pre_post_bams3 }

process samtools_stats {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

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
    mkdir -p samtools_stats/aligned samtools_stats/trimmed logs/samtools_stats
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
    mkdir -p samtools_coverage/aligned samtools_coverage/trimmed logs/samtools_coverage
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

  input:
  set val(sample), file(aligned), file(trimmed) from pre_post_bams3

  when:
  params.samtools_flagstat

  output:
  file("samtools_flagstat/{aligned,trimmed}/${sample}.flagstat.txt")
  file("logs/samtools_flagstat/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p samtools_flagstat/aligned samtools_flagstat/trimmed logs/samtools_flagstat
    log_file=logs/samtools_flagstat/!{sample}.!{workflow.sessionId}.log
    err_file=logs/samtools_flagstat/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    samtools flagstat !{aligned} 2>> $err_file > samtools_flagstat/aligned/!{sample}.flagstat.txt
    samtools flagstat !{trimmed} 2>> $err_file > samtools_flagstat/trimmed/!{sample}.flagstat.txt
  '''
}

kraken2_db = params.kraken2 ? Channel.fromPath(params.kraken2_db, type:'dir').view { "Kraken2 database : $it" } : Channel.empty()

clean_reads_classification
  .combine(kraken2_db)
  .set{ clean_reads_kraken2 }

process kraken2 {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus params.maxcpus

  when:
  params.kraken2

  input:
  set val(sample), file(clean), val(paired_single), path(kraken2_db) from clean_reads_kraken2

  output:
  file("kraken2/${sample}_kraken2_report.txt")
  file("logs/kraken2/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(percentage_cov) into kraken2_sars_results
  tuple sample, env(percentage_human) into kraken2_human_results

  shell:
  '''
    mkdir -p kraken2 logs/kraken2
    log_file=logs/kraken2/!{sample}.!{workflow.sessionId}.log
    err_file=logs/kraken2/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    kraken2 --version >> $log_file

    if [ ! -d !{kraken2_db} ]
    then
      echo "Kraken2 database could not be found. Please specify with params.kraken2_db" | tee -a $err_file
    fi

    if [ "!{paired_single}" == "single" ]
    then
      kraken2 \
        --classified-out cseqs#.fq \
        --threads !{task.cpus} \
        --db !{kraken2_db} \
        !{clean} \
        --report kraken2/!{sample}_kraken2_report.txt \
        2>> $err_file >> $log_file
    else
      kraken2 --paired \
        --classified-out cseqs#.fq \
        --threads !{task.cpus} \
        --db !{kraken2_db} \
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

process bedtools {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  when:
  params.bedtools

  input:
  set val(sample), file(bam), file(bai), file(primer_bed) from trimmed_bam_bai

  output:
  file("bedtools/${sample}.multicov.txt")
  tuple sample, env(num_failed_amplicons) into bedtools_results
  file("logs/bedtools/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p bedtools logs/bedtools
    log_file=logs/bedtools/!{sample}.!{workflow.sessionId}.log
    err_file=logs/bedtools/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    bedtools --version >> $log_file

    cat !{primer_bed} | \
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

  when:
  params.samtools_ampliconstats

  input:
  set val(sample), file(bam), file(primer_bed) from trimmed_bams_ampliconstats

  output:
  file("samtools_ampliconstats/${sample}_ampliconstats.txt")
  tuple sample, env(num_failed_amplicons) into samtools_ampliconstats_results
  file("logs/samtools_ampliconstats/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p samtools_ampliconstats logs/samtools_ampliconstats
    log_file=logs/samtools_ampliconstats/!{sample}.!{workflow.sessionId}.log
    err_file=logs/samtools_ampliconstats/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    samtools ampliconstats !{primer_bed} !{bam} 2>> $err_file > samtools_ampliconstats/!{sample}_ampliconstats.txt

    num_failed_amplicons=$(grep ^FREADS samtools_ampliconstats/!{sample}_ampliconstats.txt | cut -f 2- | tr '\t' '\n' | awk '{ if ($1 < 20) print $0 }' | wc -l)
    if [ -z "$num_failed_amplicons" ] ; then num_failed_amplicons=0 ; fi
  '''
}

process pangolin {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus params.maxcpus

  when:
  params.pangolin

  input:
  set val(sample), file(fasta) from consensus_pangolin

  output:
  file("pangolin/${sample}/lineage_report.csv")
  tuple sample, env(pangolin_lineage) into pangolin_lineage_results
  tuple sample, env(pangolin_status) into pangolin_status_results
  file("logs/pangolin/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p pangolin logs/pangolin
    log_file=logs/pangolin/!{sample}.!{workflow.sessionId}.log
    err_file=logs/pangolin/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    pangolin --version >> $log_file

    pangolin --outdir pangolin/!{sample} !{fasta} 2>> $err_file >> $log_file

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
  cpus params.medcpus

  when:
  params.nextclade

  input:
  set val(sample), file(fasta) from consensus_nextclade

  output:
  file("nextclade/${sample}_nextclade_report.csv")
  tuple sample, env(nextclade_clade) into nextclade_clade_results
  file("logs/nextclade/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p nextclade logs/nextclade
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

process summary {
  publishDir "${params.outdir}", mode: 'copy', overwrite: true
  tag "${sample}"
  echo false
  cpus 1

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
    mkdir -p summary logs/summary
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

  input:
  file(summary) from summary.collect()

  output:
  file("summary.txt")
  file("run_results.txt")
  file("logs/summary/summary.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p submission_files logs/summary
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
  cpus params.maxcpus

  input:
  file(consensus) from qc_consensus_15000_mafft.collect()
  file(reference_genome) from reference_genome_mafft

  output:
  file("mafft/mafft_aligned.fasta") into msa_file
  file("mafft/mafft_aligned.fasta") into msa_file2
  file("logs/mafft/mafft.${workflow.sessionId}.{log,err}")

  when:
  params.relatedness

  shell:
  '''
    mkdir -p mafft logs/mafft
    log_file=logs/mafft/mafft.!{workflow.sessionId}.log
    err_file=logs/mafft/mafft.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    echo "mafft version:" >> $log_file
    mafft --version 2>&1 >> $log_file

    echo ">!{params.outgroup}" > reference.fasta
    grep -v ">" !{reference_genome} >> reference.fasta

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
  tag "createing snp matrix with snp-dists"
  echo false
  cpus params.medcpus

  when:
  params.snpdists

  input:
  file(msa) from msa_file

  output:
  file("snp-dists/snp-dists.txt")
  file("logs/snp-dists/snp-dists.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p snp-dists logs/snp-dists
    log_file=logs/snp-dists/snp-dists.!{workflow.sessionId}.log
    err_file=logs/snp-dists/snp-dists.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    snp-dists -v >> $log_file

    snp-dists !{msa} > snp-dists/snp-dists.txt 2> $err_file
  '''
}

process iqtree {
  publishDir "${params.outdir}", mode: 'copy'
  tag "Creating phylogenetic tree with iqtree"
  echo false
  cpus params.maxcpus

  when:
  params.iqtree

  input:
  file(msa) from msa_file2

  output:
  file("iqtree/iqtree.{iqtree,treefile,mldist,log}")
  file("logs/iqtree/iqtree.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p iqtree logs/iqtree
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

fastq_reads_rename
   .join(consensus_rename, by:0)
   .combine(sample_file)
   .set { rename }

process rename {
  publishDir "${params.outdir}", mode: 'copy'
  tag "Renaming files for ${sample}"
  echo false
  cpus 1

  input:
  set val(sample), file(reads), val(paired_single), file(consensus), val(num_ACTG), file(sample_file) from rename

  when:
  params.sample_file.exists() && params.rename

  output:
  file("submission_files/*{genbank,gisaid}.fa") optional true into submission_fastas
  file("submission_files/*.fastq.gz")
  file("logs/rename/${sample}.${workflow.sessionId}.{err,log}")

  shell:
  '''
    mkdir -p submission_files logs/rename
    log_file=logs/rename/!{sample}.!{workflow.sessionId}.log
    err_file=logs/rename/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null

    sample_id_column=$(head -n 1 !{sample_file} | tr ',' '\\n' | grep -inw "Sample_ID" | cut -f 1 -d ':' )
    submission_id_column=$(head -n 1 !{sample_file} | tr ',' '\\n' | grep -inw "Submission_id" | cut -f 1 -d ':' )
    collection_date_column=$(head -n 1 !{sample_file} | tr ',' '\\n' | grep -inw "Collection_Date" | cut -f 1 -d ':'  )

    if [ -z $sample_id_column ] || [ -z "$submission_id_column" ] || [ -z "$collection_date_column" ]
    then
      echo "!{params.sample_file} is not the correct format"
      echo "Sorry to be overly picky, but this file needs to be a plain text file with values separated by commas (and no commas in the values)"
      echo "Required headers are 'Sample_ID','Submission_ID','Collection_Date'"
      echo "Please read documentation at https://github.com/StaPH-B/staphb_toolkit/tree/master/staphb_toolkit/workflows/cecret"
      exit 1
    fi

    all_samples=($(cut -f $sample_id_column -d "," !{sample_file}  ))
    sample_id_check=1
    for potential_sample in ${all_samples[@]}
    do
      sample_id_check=$(echo !{sample} | grep $potential_sample | head -n 1)
      if [ -n "$sample_id_check" ]
      then
        sample_id=$potential_sample
      fi
    done

    if [ -n "$sample_id" ] ; then sample_line=$(cat !{sample_file} | grep -w $sample_id | head -n 1) ; fi
    if [ -n "$sample_line" ]
    then
      echo "The line from !{params.sample_file} corresponding to !{sample} is $sample_line" | tee -a $log_file
      sample_line=$(echo $sample_line | sed 's/","/,/g')
      submission_id=$(echo $sample_line | cut -f $submission_id_column -d ',')
      if [ -z "$submission_id" ] ; then submission_id=$sample_id ; fi
      collection_date=$(echo $sample_line | cut -f $collection_date_column -d ',')
      if [ -z "$collection_date" ] ; then collection_date="missing" ; fi

      sample_file_header_reduced=$(head -n 1 !{sample_file} | tr "," '\\n' | grep -iv "Sample_ID" | grep -iv "Collection_Date" | grep -vi "Submission_ID" | tr '\\n' ' ' )
      genbank_fasta_header=">$submission_id "

      for column in ${sample_file_header_reduced[@]}
      do
        column_number=$(head -n 1 !{sample_file} | tr "," "\\n" | grep -n "$column" | cut -f 1 -d ':')
        column_value=$(echo $sample_line | cut -f $column_number -d ',')
        if [ -z "$column_value" ] ; then column_value="missing" ; fi
        genbank_fasta_header=$genbank_fasta_header"["$column"="$column_value"]"
      done

      if [ "$collection_date" == "missing" ]
      then
        year=$(date "+%Y")
        echo "The collection date is $collection_date for !{sample}" | tee -a $log_file
      else
        collection_date=$(date -d "$collection_date" "+%Y-%m-%d") || echo "Invalid date format. Try something like yyyy-mm-dd and '-resume' the workflow."
        year=$(date -d "$collection_date" "+%Y")
        echo "The collection date is $collection_date for !{sample}" | tee -a $log_file
        genbank_fasta_header=$genbank_fasta_header"[Collection_Date="$collection_date"]"
      fi

      country_check=$(echo $sample_file_header_reduced | grep -wi "country" | head -n 1 )
      if [ -z "$country_check" ]
      then
        genbank_fasta_header=$genbank_fasta_header"[Country=USA]"
        country="USA"
      else
        column_number=$(head -n 1 !{sample_file} | tr "," "\\n" | grep -in "country" | cut -f 1 -d ':')
        country=$(echo $sample_line | cut -f $column_number -d ',')
        if [ -z "$country" ] ; then country="missing" ; fi
      fi

      host_check=$(echo $sample_file_header_reduced | grep -wi "host"  | head -n 1 )
      if [ -z "$host_check" ]
      then
        genbank_fasta_header=$genbank_fasta_header"[Host=Human]"
        host="Human"
      else
        column_number=$(head -n 1 !{sample_file} | tr "," "\\n" | grep -in "host" | cut -f 1 -d ':')
        host=$(echo $sample_line | cut -f $column_number -d ',')
        if [ -z "$host" ] ; then host="missing" ; fi
      fi

      isolate_check=$(echo $sample_file_header_reduced | grep -wi "isolate"  | head -n 1 )
      if [ -z "$isolate_check" ]
      then
        organism_check=$(head -n 1 !{sample_file} | tr ',' '\\n' | grep -i "organism" | head -n 1 )
        if [ -z "$organism_check" ]
        then
          genbank_organism='SARS-CoV-2'
          gisaid_organism='hCoV-19'
        else
          column_number=$(head -n 1 !{sample_file} | tr "," "\\n" | grep -in "organism" | cut -f 1 -d ':')
          genbank_organism=$(echo $sample_line | cut -f $column_number -d ',')
          gisaid_organism=$(echo $sample_line  | cut -f $column_number -d ',')
          if [ -z "$genbank_organism" ] ; then genbank_organism="missing" ; fi
          if [ -z "$gisaid_organism" ] ; then gisaid_organism="missing" ; fi
        fi
        genbank_fasta_header=$genbank_fasta_header"[Isolate="$genbank_organism"/"$host"/"$country"/"$submission_id"/"$year"]"
      fi

      gisaid_fasta_header=">$gisaid_organism/$country/$submission_id/$year"

      if [ "!{num_ACTG}" -gt "!{params.gisaid_threshold}" ]
      then
        echo $gisaid_fasta_header > submission_files/$submission_id.gisaid.fa
        grep -v ">" !{consensus} | fold -w 75 >> submission_files/$submission_id.genbank.fa
      fi

      if [ "!{num_ACTG}" -gt "!{params.genbank_threshold}" ]
      then
        echo $genbank_fasta_header > submission_files/$submission_id.genbank.fa
        grep -v ">" !{consensus} | sed 's/^N*N//g' | fold -w 75 >> submission_files/$submission_id.genbank.fa
      fi

      if [ "!{paired_single}" == "single" ]
      then
        cp !{reads[0]} submission_files/$submission_id.fastq.gz
      else
        final_fastq_R1="submission_files/"$submission_id"_R1.fastq.gz"
        final_fastq_R2="submission_files/"$submission_id"_R2.fastq.gz"
        cp !{reads[0]} "$final_fastq_R1"
        cp !{reads[1]} "$final_fastq_R2"
      fi
    else
      echo "!{sample} was not found in !{sample_file}" | tee -a $log_file $err_file
    fi
  '''
}

process combine_fastas {
  publishDir "${params.outdir}", mode: 'copy'
  tag "Combining fastas into one multifasta"
  echo false
  cpus 1

  input:
  file(fastas) from submission_fastas.collect()

  output:
  file("submission_files/genbank_submission*.fasta") optional true
  file("submission_files/gisaid_submission*.fasta") optional true
  file("logs/combine_genbank/combine_genbank.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p submission_files logs/combine_genbank
    log_file=logs/combine_genbank/combine_genbank.!{workflow.sessionId}.log
    err_file=logs/combine_genbank/combine_genbank.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null

    cat *genbank.fa | tee submission_files/genbank_submission.fasta submission_files/genbank_submission_!{workflow.sessionId}.fasta
    cat *gisaid.fa  | tee submission_files/gisaid_submission.fasta submission_files/gisaid_submission_!{workflow.sessionId}.fasta
  '''
}

workflow.onComplete {
    println("Pipeline completed at: $workflow.complete")
    println("A summary of results can be found in a tab-delimited file: ${workflow.launchDir}/run_results.txt")
    println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
}
