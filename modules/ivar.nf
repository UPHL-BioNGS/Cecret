process ivar_consensus {
  tag           "${sample}"
  label         "process_medium"
  memory        { 2.GB * task.attempt }
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/ivar:1.4.2'

  //#UPHLICA maxForks      10
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  input:
  tuple val(sample), file(bam), file(reference_genome)

  output:
  path "consensus/${sample}.consensus.fa",                                                            emit: consensus
  path "ivar_consensus/${sample}.consensus.qual.txt",                                                 emit: qual
  path "ivar_consensus/*"
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"
  tuple val("ivar consensus"), env(ivar_version),                                                     emit: ivar_version

  shell:
  '''
    mkdir -p ivar_consensus consensus logs/!{task.process}
    log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    date > $log
    samtools --version >> $log
    ivar version >> $log
    ivar_version=$(ivar version | grep "version")

    samtools mpileup -A -d !{params.mpileup_depth} -B -Q 0 --reference !{reference_genome} !{bam} | \
      ivar consensus !{params.ivar_consensus_options} -m !{params.minimum_depth} -p ivar_consensus/!{sample}.consensus | tee -a $log

    if [ -f "ivar_consensus/!{sample}.consensus.fa" ]
    then
      echo ">!{sample}"                                               > consensus/!{sample}.consensus.fa
      grep -v ">" ivar_consensus/!{sample}.consensus.fa | fold -w 75 >> consensus/!{sample}.consensus.fa
    fi
  '''
}

process ivar_variants {
  tag           "${sample}"
  label         "process_medium"
  memory        { 2.GB * task.attempt }
  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore'}
  publishDir    "${params.outdir}", mode: 'copy'  
  container     'staphb/ivar:1.4.2'

  //#UPHLICA maxForks 10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  when:
  params.ivar_variants

  input:
  tuple val(sample), file(bam), file(reference_genome), file(gff_file)

  output:
  path "ivar_variants/${sample}.variants.tsv",      emit: variant_tsv
  path "ivar_variants/${sample}.ivar_variants.vcf", emit: ivar_variant_file
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"

  shell:
  '''
    mkdir -p ivar_variants logs/!{task.process}
    log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log
    samtools --version >> $log
    ivar version >> $log

    samtools mpileup -A -d !{params.mpileup_depth} -B -Q 0 --reference !{reference_genome} !{bam} | \
      ivar variants -p ivar_variants/!{sample}.variants !{params.ivar_variants_options} -m !{params.minimum_depth} -r !{reference_genome} -g !{gff_file} | tee -a $log

    echo '##fileformat=VCFv4.2'                                                                               >  ivar_variants/!{sample}.ivar_variants.vcf
    echo '##source=iVar'                                                                                      >> ivar_variants/!{sample}.ivar_variants.vcf
    echo '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">'                                     >> ivar_variants/!{sample}.ivar_variants.vcf
    echo '##FILTER=<ID=PASS,Description="Result of p-value <= 0.05">'                                         >> ivar_variants/!{sample}.ivar_variants.vcf
    echo '##FILTER=<ID=FAIL,Description="Result of p-value > 0.05">'                                          >> ivar_variants/!{sample}.ivar_variants.vcf
    echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'                                       >> ivar_variants/!{sample}.ivar_variants.vcf
    echo '##FORMAT=<ID=REF_DP,Number=1,Type=Integer,Description="Depth of reference base">'                   >> ivar_variants/!{sample}.ivar_variants.vcf
    echo '##FORMAT=<ID=REF_RV,Number=1,Type=Integer,Description="Depth of reference base on reverse reads">'  >> ivar_variants/!{sample}.ivar_variants.vcf
    echo '##FORMAT=<ID=REF_QUAL,Number=1,Type=Integer,Description="Mean quality of reference base">'          >> ivar_variants/!{sample}.ivar_variants.vcf
    echo '##FORMAT=<ID=ALT_DP,Number=1,Type=Integer,Description="Depth of alternate base">'                   >> ivar_variants/!{sample}.ivar_variants.vcf
    echo '##FORMAT=<ID=ALT_RV,Number=1,Type=Integer,Description="Deapth of alternate base on reverse reads">' >> ivar_variants/!{sample}.ivar_variants.vcf
    echo '##FORMAT=<ID=ALT_QUAL,Number=1,Type=String,Description="Mean quality of alternate base">'           >> ivar_variants/!{sample}.ivar_variants.vcf
    echo '##FORMAT=<ID=ALT_FREQ,Number=1,Type=String,Description="Frequency of alternate base">'              >> ivar_variants/!{sample}.ivar_variants.vcf
    echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t!{sample}'                                >> ivar_variants/!{sample}.ivar_variants.vcf
    tail -n+2 !{task.process}/!{sample}.variants.tsv | \
      awk '{print $1 "\t" $2 "\t.\t" $3 "\t" $4 "\t.\t.\tREF_DP=" $5 ";REF_RV=" $6 ";REF_QUAL=" $7 ";ALT_DP=" $8 ";ALT_RV=" $9 ";ALT_QUAL=" $10 "\tGT:PL\t1/1:" $12 "," $12-$8 "," $8 }' \
      >> ivar_variants/!{sample}.ivar_variants.vcf
  '''
}

process ivar_trim {
  tag        "${sample}"
  label      "process_medium"
  publishDir "${params.outdir}", mode: 'copy'
  container  'staphb/ivar:1.4.2'
  
  //#UPHLICA maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '45m'

  input:
  tuple val(sample), file(bam), file(primer_bed)

  output:
  tuple val(sample), file("ivar_trim/${sample}.primertrim.sorted.bam"),                                                         emit: trimmed_bam
  tuple val(sample), file("ivar_trim/${sample}.primertrim.sorted.bam"), file("ivar_trim/${sample}.primertrim.sorted.bam.bai"),  emit: bam_bai
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"
  path "ivar_trim/${sample}_ivar.log",                                                                                          emit: ivar_trim_files
  tuple val("${params.trimmer}"), env(trimmer_version),                                                                         emit: trimmer_version

  shell:
  '''
    mkdir -p ivar_trim logs/!{task.process}
    log=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log
    ivar version >> $log
    trimmer_version="ivar : $(ivar version | grep version)"
    
    # trimming the reads
    ivar trim !{params.ivar_trim_options} -e -i !{bam} -b !{primer_bed} -p ivar_trim/!{sample}.primertrim | tee -a $log

    if [ -S "$log" ] ; then grep "Found" -A 10000 $log | grep -A 10000 "primers in BED file" 2> ivar_trim/!{sample}_ivar.log ; else touch ivar_trim/!{sample}_ivar.log ; fi
    
    # sorting and indexing the trimmed bams
    samtools sort ivar_trim/!{sample}.primertrim.bam -o ivar_trim/!{sample}.primertrim.sorted.bam | tee -a $log
    samtools index ivar_trim/!{sample}.primertrim.sorted.bam | tee -a $log
  '''
}
