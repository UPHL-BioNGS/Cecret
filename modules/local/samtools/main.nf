process SAMTOOLS_QC {
  tag        "${meta.id}"
  label      "process_single"
  container  'staphb/samtools:1.22.1'

  input:
  tuple val(meta), file(bam)

  output:
  val(meta), emit: meta // for linter
  path "samtools/*.stats.txt",       emit: stats
  path "samtools/*",                 emit: files
  path "samtools/*.cov.txt",         emit: coverage
  path "samtools/*.flagstat.txt",    emit: flagstat
  path "samtools/*.depth.txt",       emit: depth
  path "versions.yml",               emit: versions
  path "logs/${task.process}/*.log", emit: log

  when:
  task.ext.when == null || task.ext.when

  script:
  def stats_args    = task.ext.stats_args     ?: "${params.samtools_stats_options}"
  def coverage_args = task.ext.coverage_args  ?: "${params.samtools_coverage_options}"
  def flagstat_args = task.ext.flagstat_args  ?: "${params.samtools_flagstat_options}"
  def depth_args    = task.ext.depth_args     ?: "${params.samtools_depth_options}"
  def prefix        = task.ext.prefix         ?: "${meta.id}"
  def which_stage   = bam.baseName.contains('primertrim') ? 'final' : 'initial'
  """
    mkdir -p samtools logs/${task.process}
    log=logs/${task.process}/${prefix}.${which_stage}.${workflow.sessionId}.log

    date > \$log
    samtools --version >> \$log

    samtools stats \
      ${stats_args} \
      ${bam} \
      > samtools/${prefix}.${which_stage}.stats.txt

    samtools coverage \
      ${coverage_args} \
      ${bam} \
      -m \
      -o \
      samtools/${prefix}.${which_stage}.cov.hist \
      | tee -a \$log
    
    samtools coverage \
      ${coverage_args} \
      ${bam} \
      | awk -v sample=${prefix} '{print sample "\\t" \$0 }' \
      | sed '0,/${prefix}/s//sample/' \
      > samtools/${prefix}.${which_stage}.cov.txt \
      | tee -a \$log

    samtools flagstat \
      ${flagstat_args} \
      ${bam} \
      | tee samtools/${prefix}.${which_stage}.flagstat.txt

    samtools depth \
      ${depth_args} \
      ${bam} \
      > samtools/${prefix}.${which_stage}.depth.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      samtools: \$(samtools --version | head -n 1 | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}

process SAMTOOLS_AMPLICONSTATS {
  tag        "${meta.id}"
  label      "process_single"
  container  'staphb/samtools:1.22.1'
  errorStrategy 'ignore'
  
  input:
  tuple val(meta), file(bam), file(primer_bed)

  output:
  tuple val(meta), file("samtools/*_ampliconstats.txt"), emit: ampliconstats
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: "${params.samtools_ampliconstats_options}"
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p samtools logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    date > \$log
    samtools --version >> \$log

    samtools ampliconstats \
      ${args} \
      ${primer_bed} \
      ${bam} > samtools/${prefix}_ampliconstats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      samtools: \$(samtools --version | head -n 1 | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}

process SAMTOOLS_PLOT_AMPLICONSTATS {
  tag           "${meta.id}"
  label         "process_single"
  container     'staphb/samtools:1.22.1'
  // fails for empty samples
  errorStrategy 'ignore'

  input:
  tuple val(meta), file(ampliconstats)

  output:
  val(meta), emit: meta // for linter
  path "samtools_plot_ampliconstats/*", emit: files, optional: true
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: "${params.samtools_plot_ampliconstats_options}"
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p samtools_plot_ampliconstats/${prefix} logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    date > \$log
    samtools --version >> \$log

    plot-ampliconstats ${args} \
      samtools_plot_ampliconstats/${prefix} \
      ${ampliconstats}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      samtools: \$(samtools --version | head -n 1 | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}

process SAMTOOLS_SORT {
  tag        "${meta.id}"
  label      "process_high"
  container  'staphb/samtools:1.22.1'

  input:
  tuple val(meta), file(sam)

  output:
  tuple val(meta), file("aligned/*.sorted.bam"), file("aligned/*.sorted.bam.bai"), emit: bam_bai
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p aligned logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log
    samtools --version >> \$log

    samtools sort -@ ${task.cpus} ${sam} | \
      samtools view -F 4 -o aligned/${prefix}.sorted.bam | tee -a \$log

    # indexing the bams
    samtools index aligned/${prefix}.sorted.bam | tee -a \$log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      samtools: \$(samtools --version | head -n 1 | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}

process SAMTOOLS_FILTER {
  tag        "${meta.id}"
  label      "process_single"
  container  'staphb/samtools:1.22.1'

  input:
  tuple val(meta), file(sam)

  output:
  tuple val(meta), file("filtered/*_filtered_{R1,R2}.fastq.gz"), optional: true, emit: filtered_reads
  path "filtered/*_filtered_unpaired.fastq.gz", optional: true, emit: unpaired
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: '-F 4'
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p filtered logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log
    samtools --version >> \$log

    samtools sort -n ${sam} | \
      samtools fastq ${args} ${params.filter_options} \
      -s filtered/${prefix}_filtered_unpaired.fastq.gz \
      -1 filtered/${prefix}_filtered_R1.fastq.gz \
      -2 filtered/${prefix}_filtered_R2.fastq.gz \
      | tee -a \$log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      samtools: \$(samtools --version | head -n 1 | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}

process SAMTOOLS_AMPLICONCLIP {
  tag        "${meta.id}"
  label      "process_single"
  container  'staphb/samtools:1.22.1'


  input:
  tuple val(meta), file(bam), file(primer_bed)

  output:
  tuple val(meta), file("ampliconclip/*.primertrim.sorted.bam"), file("ampliconclip/*.primertrim.sorted.bam.bai"), emit: bam_bai
  path "logs/${task.process}/*.log", emit: log                                                         
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: "${params.samtools_ampliconclip_options}"
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p ampliconclip logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log
    samtools --version >> \$log
    trimmer_version=\$(samtools --version | head -n 1 | awk '{print \$NF}')

    # trimming the reads
    samtools ampliconclip ${args} -b ${primer_bed} ${bam} | \
      samtools sort |  \
      samtools view -F 4 -o ampliconclip/${prefix}.primertrim.sorted.bam | tee -a \$log

    samtools index ampliconclip/${prefix}.primertrim.sorted.bam | tee -a \$log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      samtools: \$(samtools --version | head -n 1 | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}

process SAMTOOLS_MARKDUP {
  tag        "${meta.id}"
  label      "process_single"
  container  'staphb/samtools:1.22.1'

  input:
  tuple val(meta), val(type), file(sam) 

  output:
  tuple val(meta), file("markdup/*.markdup.sorted.bam"), file("markdup/*.markdup.sorted.bam.bai"), emit: bam_bai
  path "markdup/*_markdupstats.txt", emit: stats
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions
  
  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: "${params.samtools_markdup_options}"
  def fmargs = task.ext.fmargs ?: "${params.samtools_fixmate_options}"
  def prefix = task.ext.prefix ?: "${meta.id}"
  if ( meta.single_end ) {
  """
    mkdir -p markdup logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log
    samtools --version >> \$log

    samtools sort ${sam} | \
      samtools markdup ${args} -@ ${task.cpus} -s -f markdup/${prefix}_markdupstats.txt - markdup/${prefix}.markdup.sorted.bam | tee -a \$log

    samtools index markdup/${prefix}.markdup.sorted.bam | tee -a \$log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      samtools: \$(samtools --version | head -n 1 | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
  } else {
  """
    mkdir -p markdup logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > \$log
    samtools --version >> \$log

    samtools sort -n ${sam} | \
      samtools fixmate ${fmargs} -m -@ ${task.cpus} - - | \
      samtools sort | \
      samtools markdup ${args} -@ ${task.cpus} -s -f markdup/${prefix}_markdupstats.txt - markdup/${prefix}.markdup.sorted.bam | tee -a \$log

    samtools index markdup/${prefix}.markdup.sorted.bam | tee -a \$log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      samtools: \$(samtools --version | head -n 1 | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
  }
}

