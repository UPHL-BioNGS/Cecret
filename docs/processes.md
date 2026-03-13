# Processes

This page documents all processes in the pipeline.

## Contents

- [ACI](#aci)
- [ARTIC](#artic)
- [BBNORM](#bbnorm)
- [BCFTOOLS](#bcftools)
- [BWA](#bwa)
- [DATASETS](#datasets)
- [ENA](#ena)
- [FASTP](#fastp)
- [FASTQC](#fastqc)
- [FREYJA](#freyja)
- [FREYJA_AGGREGATE](#freyja-aggregate)
- [FREYJA_PATHOGEN](#freyja-pathogen)
- [FREYJA_UPDATE](#freyja-update)
- [HEATCLUSTER](#heatcluster)
- [IGV_REPORTS](#igv-reports)
- [IQTREE](#iqtree)
- [IVAR_CONSENSUS](#ivar-consensus)
- [IVAR_VARIANTS](#ivar-variants)
- [IVAR_TRIM](#ivar-trim)
- [KRAKEN2](#kraken2)
- [KRAKEN2_DB](#kraken2-db)
- [PREP](#prep)
- [SUMMARY](#summary)
- [UNZIP](#unzip)
- [MAFFT](#mafft)
- [MINIMAP2](#minimap2)
- [MULTIQC](#multiqc)
- [NEXTCLADE_DATASET](#nextclade-dataset)
- [NEXTCLADE](#nextclade)
- [PANGO_ALIASOR](#pango-aliasor)
- [PANGOLIN](#pangolin)
- [PHYTREEVIZ](#phytreeviz)
- [SAMTOOLS_QC](#samtools-qc)
- [SAMTOOLS_AMPLICONSTATS](#samtools-ampliconstats)
- [SAMTOOLS_PLOT_AMPLICONSTATS](#samtools-plot-ampliconstats)
- [SAMTOOLS_SORT](#samtools-sort)
- [SAMTOOLS_FILTER](#samtools-filter)
- [SAMTOOLS_AMPLICONCLIP](#samtools-ampliconclip)
- [SAMTOOLS_MARKDUP](#samtools-markdup)
- [SEQYCLEAN](#seqyclean)
- [SNPDISTS](#snpdists)
- [VADR](#vadr)
- [ARTIC_FILTER](#artic-filter)

## ACI {#aci}

*Defined in `modules/local/aci/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(bam), file(bed)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `meta` | `val` | `meta` | - |
| `aci/*/*_amplicon_depth.csv` | `path` | `cov` | - |
| `aci/*/*_amplicon_depth.png` | `path` | `for_multiqc` | - |
| `aci/*/*` | `path` | `everything` | - |
| `logs/${task.process` | `path` | - | - |


## ARTIC {#artic}

*Defined in `modules/local/artic/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(fastq), file(reference), file(bed)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta), file("artic/*.primertrim.sorted.bam"), file("artic/*.primertrim.sorted.bam.bai")` | `tuple` | `bam` | - |
| `consensus/*.consensus.fa` | `path` | `consensus` | - |
| `artic/*vcf*` | `path` | `vcf` | - |
| `artic/*{tsv` | `path` | - | - |


## BBNORM {#bbnorm}

*Defined in `modules/local/bbnorm/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(reads)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta), path("bbnorm/*.fastq.gz")` | `tuple` | `fastq` | - |
| `logs/*/*.log` | `path` | `log` | - |
| `versions.yml` | `path` | `versions` | - |


## BCFTOOLS {#bcftools}

*Defined in `modules/local/bcftools/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(bam), file(reference_genome)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta), file("bcftools_variants/*.vcf")` | `tuple` | `vcf` | - |
| `bcftools_variants/*.vcf` | `path` | `bcftools_variants_file` | - |
| `logs/${task.process` | `path` | - | - |


## BWA {#bwa}

*Defined in `modules/local/bwa/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(reads), file(reference_genome)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta), file("bwa/*.sam")` | `tuple` | `sam` | - |
| `logs/${task.process` | `path` | - | - |


## DATASETS {#datasets}

*Defined in `modules/local/datasets/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `genomes/*fasta` | `path` | `fasta` | - |
| `versions.yml` | `path` | `versions` | - |


## ENA {#ena}

*Defined in `modules/local/ena/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(SRR), file("ena/paired/*fastq.gz"), val(false)` | `tuple` | `paired` | - |
| `val(SRR), file("ena/single/*fastq.gz"), val(true)` | `tuple` | `single` | - |
| `logs/*/*.log` | `path` | `log` | - |
| `versions.yml` | `path` | `versions` | - |


## FASTP {#fastp}

*Defined in `modules/local/fastp/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(reads)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta)` | `tuple` | - | - |


## FASTQC {#fastqc}

*Defined in `modules/local/fastqc/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(fastq)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `meta` | `val` | `meta` | - |
| `fastqc/*.html` | `path` | `files` | - |
| `fastqc/*_fastqc.zip` | `path` | `fastqc_files` | - |
| `fastqc/*_fastq_name.csv` | `path` | `fastq_name` | - |
| `logs/${task.process` | `path` | - | - |


## FREYJA {#freyja}

*Defined in `modules/local/freyja/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(bam), file(reference_genome)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta)` | `tuple` | - | - |


## FREYJA_AGGREGATE {#freyja-aggregate}

*Defined in `modules/local/freyja/main.nf:51`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `freyja/*` | `path` | `files` | - |
| `freyja/aggregated-freyja.tsv` | `path` | `aggregated_freyja_file` | - |
| `freyja/*mqc.png` | `path` | `for_multiqc` | - |
| `logs/${task.process` | `path` | - | - |


## FREYJA_PATHOGEN {#freyja-pathogen}

*Defined in `modules/local/freyja/main.nf:106`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(bam), file(reference_genome), file(db)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta)` | `tuple` | - | - |


## FREYJA_UPDATE {#freyja-update}

*Defined in `modules/local/freyja/main.nf:159`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `db` | `path` | `db` | - |
| `logs/${task.process` | `path` | - | - |


## HEATCLUSTER {#heatcluster}

*Defined in `modules/local/heatcluster/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `heatcluster/*` | `path` | `files` | - |
| `heatcluster/*.png` | `path` | `for_multiqc` | - |
| `logs/${task.process` | `path` | - | - |


## IGV_REPORTS {#igv-reports}

*Defined in `modules/local/igvreports/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(vcf), file(bam), file(bai), file(reference_genome)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `meta` | `val` | `meta` | - |
| `igv_reports/*` | `path` | `report` | - |
| `logs/${task.process` | `path` | - | - |


## IQTREE {#iqtree}

*Defined in `modules/local/iqtree/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `iqtree/*.{iqtree` | `path` | - | - |


## IVAR_CONSENSUS {#ivar-consensus}

*Defined in `modules/local/ivar/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(bam), file(reference_genome)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `meta` | `val` | `meta` | - |
| `consensus/*.consensus.fa` | `path` | `consensus` | - |
| `ivar_consensus/*.consensus.qual.txt` | `path` | `qual` | - |
| `ivar_consensus/*` | `path` | `everything` | - |
| `logs/${task.process` | `path` | - | - |


## IVAR_VARIANTS {#ivar-variants}

*Defined in `modules/local/ivar/main.nf:49`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(bam), file(reference_genome), file(gff_file)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `meta` | `val` | `meta` | - |
| `ivar_variants/*.variants.tsv` | `path` | `variant_tsv` | - |
| `ivar_variants/*.ivar_variants.vcf` | `path` | `vcf` | - |
| `logs/${task.process` | `path` | - | - |


## IVAR_TRIM {#ivar-trim}

*Defined in `modules/local/ivar/main.nf:108`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(bam), file(primer_bed)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta), file("ivar_trim/*.primertrim.sorted.bam")` | `tuple` | `trimmed_bam` | - |
| `val(meta), file("ivar_trim/*.primertrim.sorted.bam"), file("ivar_trim/*.primertrim.sorted.bam.bai")` | `tuple` | `bam_bai` | - |
| `logs/${task.process` | `path` | - | - |


## KRAKEN2 {#kraken2}

*Defined in `modules/local/kraken2/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(clean)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `meta` | `val` | `meta` | - |
| `kraken2/*_kraken2_report.txt` | `path` | `kraken2_files` | - |
| `kraken2/*` | `path` | `everything` | - |
| `logs/${task.process` | `path` | - | - |


## KRAKEN2_DB {#kraken2-db}

*Defined in `modules/local/kraken2/main.nf:48`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(clean), path(kraken2_db)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `meta` | `val` | `meta` | - |
| `kraken2/*_kraken2_report.txt` | `path` | `kraken2_files` | - |
| `kraken2/*` | `path` | `everything` | - |
| `logs/${task.process` | `path` | - | - |


## PREP {#prep}

*Defined in `modules/local/local/main.nf:2`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(fasta)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `fasta_prep/*` | `path` | `fastas` | - |
| `versions.yml` | `path` | `versions` | - |


## SUMMARY {#summary}

*Defined in `modules/local/local/main.nf:33`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `file(files), file(script), file(multiqc)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `cecret_results.{csv` | `path` | - | - |


## UNZIP {#unzip}

*Defined in `modules/local/local/main.nf:69`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `dataset` | `path` | `dataset` | - |
| `versions.yml` | `path` | `versions` | - |


## MAFFT {#mafft}

*Defined in `modules/local/mafft/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `mafft/mafft_aligned.fasta` | `path` | `msa` | - |
| `mafft/*` | `path` | `files` | - |
| `logs/${task.process` | `path` | - | - |


## MINIMAP2 {#minimap2}

*Defined in `modules/local/minimap2/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(reads), file(reference_genome)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta), file("aligned/*.sam")` | `tuple` | `sam` | - |
| `logs/${task.process` | `path` | - | - |


## MULTIQC {#multiqc}

*Defined in `modules/local/multiqc/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `multiqc/multiqc_report.html` | `path` | `html` | - |
| `multiqc/multiqc_data/*` | `path` | `files` | - |
| `multiqc/multiqc_data` | `path` | `multiqc_data` | - |
| `software_versions.yml` | `path` | `versions` | - |
| `logs/${task.process` | `path` | - | - |


## NEXTCLADE_DATASET {#nextclade-dataset}

*Defined in `modules/local/nextclade/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `dataset` | `path` | `dataset` | - |
| `logs/${task.process` | `path` | - | - |


## NEXTCLADE {#nextclade}

*Defined in `modules/local/nextclade/main.nf:38`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `nextclade/nextclade.csv` | `path` | `nextclade_file` | - |
| `nextclade/*` | `path` | `results` | - |
| `file("nextclade/nextclade.aligned.fasta"), file("nextclade/nextclade.nwk")` | `tuple` | `prealigned` | - |
| `logs/${task.process` | `path` | - | - |


## PANGO_ALIASOR {#pango-aliasor}

*Defined in `modules/local/pango_aliasor/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `pango_aliasor/*.csv` | `path` | `results` | - |
| `logs/${task.process` | `path` | - | - |


## PANGOLIN {#pangolin}

*Defined in `modules/local/pangolin/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `pangolin/*` | `path` | `results` | - |
| `pangolin/lineage_report.csv` | `path` | `pangolin_file` | - |
| `logs/${task.process` | `path` | - | - |


## PHYTREEVIZ {#phytreeviz}

*Defined in `modules/local/phytreeviz/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `phytreeviz/tree.png` | `path` | `for_multiqc` | - |
| `logs/${task.process` | `path` | - | - |


## SAMTOOLS_QC {#samtools-qc}

*Defined in `modules/local/samtools/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(bam)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `meta` | `val` | `meta` | - |
| `samtools/*.stats.txt` | `path` | `stats` | - |
| `samtools/*` | `path` | `files` | - |
| `samtools/*.cov.txt` | `path` | `coverage` | - |
| `samtools/*.flagstat.txt` | `path` | `flagstat` | - |
| `samtools/*.depth.txt` | `path` | `depth` | - |
| `versions.yml` | `path` | `versions` | - |
| `logs/${task.process` | `path` | - | - |


## SAMTOOLS_AMPLICONSTATS {#samtools-ampliconstats}

*Defined in `modules/local/samtools/main.nf:75`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(bam), file(primer_bed)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta), file("samtools/*_ampliconstats.txt")` | `tuple` | `ampliconstats` | - |
| `logs/${task.process` | `path` | - | - |


## SAMTOOLS_PLOT_AMPLICONSTATS {#samtools-plot-ampliconstats}

*Defined in `modules/local/samtools/main.nf:115`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(ampliconstats)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `meta` | `val` | `meta` | - |
| `samtools_plot_ampliconstats/*` | `path` | `files` | - |
| `logs/${task.process` | `path` | - | - |


## SAMTOOLS_SORT {#samtools-sort}

*Defined in `modules/local/samtools/main.nf:156`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(sam)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta), file("aligned/*.sorted.bam"), file("aligned/*.sorted.bam.bai")` | `tuple` | `bam_bai` | - |
| `logs/${task.process` | `path` | - | - |


## SAMTOOLS_FILTER {#samtools-filter}

*Defined in `modules/local/samtools/main.nf:196`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(sam)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta)` | `tuple` | - | - |


## SAMTOOLS_AMPLICONCLIP {#samtools-ampliconclip}

*Defined in `modules/local/samtools/main.nf:239`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(bam), file(primer_bed)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta), file("ampliconclip/*.primertrim.sorted.bam"), file("ampliconclip/*.primertrim.sorted.bam.bai")` | `tuple` | `bam_bai` | - |
| `logs/${task.process` | `path` | - | - |


## SAMTOOLS_MARKDUP {#samtools-markdup}

*Defined in `modules/local/samtools/main.nf:283`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), val(type), file(sam)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta), file("markdup/*.markdup.sorted.bam"), file("markdup/*.markdup.sorted.bam.bai")` | `tuple` | `bam_bai` | - |
| `markdup/*_markdupstats.txt` | `path` | `stats` | - |
| `logs/${task.process` | `path` | - | - |


## SEQYCLEAN {#seqyclean}

*Defined in `modules/local/seqyclean/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(reads)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta)` | `tuple` | - | - |


## SNPDISTS {#snpdists}

*Defined in `modules/local/snp-dists/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `snp-dists/*.txt` | `path` | `matrix` | - |
| `logs/${task.process` | `path` | - | - |


## VADR {#vadr}

*Defined in `modules/local/vadr/main.nf:1`*

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `vadr/*` | `path` | `vadr_files` | - |
| `vadr/vadr.vadr.sqa` | `path` | `vadr_file` | - |
| `logs/${task.process` | `path` | - | - |


## ARTIC_FILTER {#artic-filter}

*Defined in `modules/local/artic_filter/main.nf:1`*

### Inputs

| Name | Type | Description |
|------|------|-------------|
| `val(meta), file(fastq)` | `tuple` | - |

### Outputs

| Name | Type | Emit | Description |
|------|------|------|-------------|
| `val(meta), file("artic/*_filtered.fastq")` | `tuple` | `fastq` | - |
| `logs/${task.process` | `path` | - | - |


---

*This pipeline was built with [Nextflow](https://nextflow.io).
Documentation generated by [nf-docs](https://github.com/ewels/nf-docs) v0.2.0 on 2026-03-12 23:00:43 UTC.*
