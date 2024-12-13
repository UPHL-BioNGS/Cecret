process IVAR_CONSENSUS {
  tag           "${meta.id}"
  label         "process_medium"
  container     'staphb/ivar:1.4.3'

  input:
  tuple val(meta), file(bam), file(reference_genome)

  output:
  val(meta), emit: meta // for linter
  path "consensus/*.consensus.fa", emit: consensus
  path "ivar_consensus/*.consensus.qual.txt", emit: qual
  path "ivar_consensus/*", emit: everything
  path "logs/${task.process}/*.log", emit: log
  tuple val("ivar consensus"), env("ivar_version"), emit: ivar_version
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: "${params.ivar_consensus_options}"
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p ivar_consensus consensus logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    date > \$log
    samtools --version >> \$log
    ivar version >> \$log
    ivar_version=\$(ivar version | head -n 1 | awk '{print \$NF}')

    samtools mpileup -A -d ${params.mpileup_depth} -B -Q 0 --reference ${reference_genome} ${bam} | \
      ivar consensus ${args} -m ${params.minimum_depth} -p ivar_consensus/${prefix}.consensus | tee -a \$log

    if [ -f "ivar_consensus/${prefix}.consensus.fa" ]
    then
      echo ">${prefix}"                                               > consensus/${prefix}.consensus.fa
      grep -v ">" ivar_consensus/${prefix}.consensus.fa | fold -w 75 >> consensus/${prefix}.consensus.fa
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      ivar: \$(ivar version | head -n 1 | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}

process IVAR_VARIANTS {
  tag           "${meta.id}"
  label         "process_medium"
  container     'staphb/ivar:1.4.3'

  input:
  tuple val(meta), file(bam), file(reference_genome), file(gff_file)

  output:
  val(meta), emit: meta // for linter
  path "ivar_variants/*.variants.tsv", emit: variant_tsv
  path "ivar_variants/*.ivar_variants.vcf", emit: vcf
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: "${params.ivar_variants_options}"
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p ivar_variants logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log
    samtools --version >> \$log
    ivar version >> \$log

    samtools mpileup -A -d ${params.mpileup_depth} -B -Q 0 --reference ${reference_genome} ${bam} | \
      ivar variants -p ivar_variants/${prefix}.variants ${args} -m ${params.minimum_depth} -r ${reference_genome} -g ${gff_file} | tee -a \$log

    echo '##fileformat=VCFv4.2'                                                                               >  ivar_variants/${prefix}.ivar_variants.vcf
    echo '##source=iVar'                                                                                      >> ivar_variants/${prefix}.ivar_variants.vcf
    echo '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">'                                     >> ivar_variants/${prefix}.ivar_variants.vcf
    echo '##FILTER=<ID=PASS,Description="Result of p-value <= 0.05">'                                         >> ivar_variants/${prefix}.ivar_variants.vcf
    echo '##FILTER=<ID=FAIL,Description="Result of p-value > 0.05">'                                          >> ivar_variants/${prefix}.ivar_variants.vcf
    echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'                                       >> ivar_variants/${prefix}.ivar_variants.vcf
    echo '##FORMAT=<ID=REF_DP,Number=1,Type=Integer,Description="Depth of reference base">'                   >> ivar_variants/${prefix}.ivar_variants.vcf
    echo '##FORMAT=<ID=REF_RV,Number=1,Type=Integer,Description="Depth of reference base on reverse reads">'  >> ivar_variants/${prefix}.ivar_variants.vcf
    echo '##FORMAT=<ID=REF_QUAL,Number=1,Type=Integer,Description="Mean quality of reference base">'          >> ivar_variants/${prefix}.ivar_variants.vcf
    echo '##FORMAT=<ID=ALT_DP,Number=1,Type=Integer,Description="Depth of alternate base">'                   >> ivar_variants/${prefix}.ivar_variants.vcf
    echo '##FORMAT=<ID=ALT_RV,Number=1,Type=Integer,Description="Deapth of alternate base on reverse reads">' >> ivar_variants/${prefix}.ivar_variants.vcf
    echo '##FORMAT=<ID=ALT_QUAL,Number=1,Type=String,Description="Mean quality of alternate base">'           >> ivar_variants/${prefix}.ivar_variants.vcf
    echo '##FORMAT=<ID=ALT_FREQ,Number=1,Type=String,Description="Frequency of alternate base">'              >> ivar_variants/${prefix}.ivar_variants.vcf
    echo -e '#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t${prefix}'                       >> ivar_variants/${prefix}.ivar_variants.vcf
    tail -n+2 ${task.process}/${prefix}.variants.tsv | \
      awk '{print \$1 "\\t" \$2 "\\t.\\t" \$3 "\\t" \$4 "\\t.\\t.\\tREF_DP=" \$5 ";REF_RV=" \$6 ";REF_QUAL=" \$7 ";ALT_DP=" \$8 ";ALT_RV=" \$9 ";ALT_QUAL=" \$10 "\\tGT:PL\\t1/1:" \$12 "," \$12-\$8 "," \$8 }' \
      >> ivar_variants/${prefix}.ivar_variants.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      ivar: \$(ivar version | head -n 1 | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}

process IVAR_TRIM {
  tag        "${meta.id}"
  label      "process_medium"
  container  'staphb/ivar:1.4.3'

  input:
  tuple val(meta), file(bam), file(primer_bed)

  output:
  tuple val(meta), file("ivar_trim/*.primertrim.sorted.bam"), emit: trimmed_bam
  tuple val(meta), file("ivar_trim/*.primertrim.sorted.bam"), file("ivar_trim/*.primertrim.sorted.bam.bai"), emit: bam_bai
  path "logs/${task.process}/*.log", emit: log
  path "ivar_trim/*_ivar.log", emit: ivar_trim_files
  tuple val("${params.trimmer}"), env("trimmer_version"), emit: trimmer_version
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: "${params.ivar_trim_options}"
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p ivar_trim logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log
    ivar version >> \$log
    trimmer_version=\$(ivar version | head -n 1 | awk '{print \$NF}')
    
    # trimming the reads
    ivar trim ${args} -e -i ${bam} -b ${primer_bed} -p ivar_trim/${prefix}.primertrim | tee -a \$log

    if [ -S "\$log" ] ; then grep "Found" -A 10000 \$log | grep -A 10000 "primers in BED file" 2> ivar_trim/${prefix}_ivar.log ; else touch ivar_trim/${prefix}_ivar.log ; fi
    
    # sorting and indexing the trimmed bams
    samtools sort ivar_trim/${prefix}.primertrim.bam -o ivar_trim/${prefix}.primertrim.sorted.bam | tee -a \$log
    samtools index ivar_trim/${prefix}.primertrim.sorted.bam | tee -a \$log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      ivar: \$(ivar version | head -n 1 | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}
