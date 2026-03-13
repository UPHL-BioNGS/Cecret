# Pipeline Inputs

This page documents all input parameters for the pipeline.

## Input/output options

### `--input` {#input}

**Type:** `string` | *Optional* | **Format:** `file-path`

same as sample_sheet


### `--sample_sheet` {#sample-sheet}

**Type:** `string` | *Optional*

sample sheet with sample, fastq_1, and fastq_2 columns


### `--reads` {#reads}

**Type:** `string` | *Optional*

Directory with paired-end fastq files


### `--single_reads` {#single-reads}

**Type:** `string` | *Optional*

Directory with single-end illumina fastq files


### `--nanopore` {#nanopore}

**Type:** `string` | *Optional*

Directory with nanopore fastq files


### `--fastas` {#fastas}

**Type:** `string` | *Optional*

input channel for fastas


### `--multifastas` {#multifastas}

**Type:** `string` | *Optional*

input channel for multifasta files


### `--sra_accessions` {#sra-accessions}

**Type:** `string` | *Optional*

list encased in [] brackets of fastq accessions to download from SRA/ENA (Illumina reads only)

**Default:** `[]`


### `--genome_accessions` {#genome-accessions}

**Type:** `string` | *Optional*

list encased in [] brackets of genome accessions to download with datasets

**Default:** `[]`


### `--species` {#species}

**Type:** `string` | *Optional*

specifies species-specific sub-workflows

**Default:** `sarscov2`

**Allowed values:**
- `other`
- `mpx`
- `sarscov2`


### `--vadr_reference` {#vadr-reference}

**Type:** `string` | *Optional*

Specifies reference for vadr in container

**Default:** `sarscov2`


### `--freyja_pathogen` {#freyja-pathogen}

**Type:** `string` | *Optional*

Specifies freyja pathogen


### `--freyja_update` {#freyja-update}

**Type:** `boolean` | *Optional*

Whether to update freyja db

**Default:** `True`


### `--nextclade_dataset` {#nextclade-dataset}

**Type:** `string` | *Optional*

Specifies nextclade dataset

**Default:** `sars-cov-2`


### `--iqtree_outgroup` {#iqtree-outgroup}

**Type:** `string` | *Optional*

outgroup for multiple sequence alignment

**Default:** `MN908947`


### `--minimum_depth` {#minimum-depth}

**Type:** `integer` | *Optional*

minimum depth for calling a variant

**Default:** `100`


### `--mpileup_depth` {#mpileup-depth}

**Type:** `integer` | *Optional*

number of reads put into memory by samtools/bcftools

**Default:** `8000`


### `--kraken2_db` {#kraken2-db}

**Type:** `string` | *Optional*

directory to kraken2 database


### `--outdir` {#outdir}

**Type:** `string` | **Required** | **Format:** `directory-path`

The output directory where the results will be saved. Absolute paths are required on cloud infrastructure.

**Default:** `cecret`



## Reference files

### `--reference_genome` {#reference-genome}

**Type:** `string` | *Optional*

THE Reference genome


### `--amplicon_bed` {#amplicon-bed}

**Type:** `string` | *Optional*

Bedfile for amplicons


### `--gff` {#gff}

**Type:** `string` | *Optional*

File used in ivar variants. Must correspond with reference genome.


### `--primer_bed` {#primer-bed}

**Type:** `string` | *Optional*

File with bedfile of primers used in the analysis


### `--primer_set` {#primer-set}

**Type:** `string` | *Optional*

Specifies a primer set included in repo

**Default:** `ncov_V5.3.2`

**Allowed values:**
- `midnight_idt_V1`
- `midnight_ont_V1`
- `midnight_ont_V2`
- `midnight_ont_V3`
- `ncov_V3`
- `ncov_V4`
- `ncov_V4.1`
- `ncov_V5.3.2`
- `mpx_primalseq`
- `mpx_idt`



## Workflow Components

### `--cleaner` {#cleaner}

**Type:** `string` | *Optional*

Specifies what tool to use to remove low quality reads

**Default:** `seqyclean`

**Allowed values:**
- `seqyclean`
- `fastp`


### `--aligner` {#aligner}

**Type:** `string` | *Optional*

Specifies which aligner is going to be used.

**Default:** `bwa`

**Allowed values:**
- `bwa`
- `minimap2`


### `--trimmer` {#trimmer}

**Type:** `string` | *Optional*

Specifies which tool to use to trim primers from primerbedfile

**Default:** `ivar`

**Allowed values:**
- `samtools`
- `ivar`
- `none`


### `--msa` {#msa}

**Type:** `string` | *Optional*

Specifies what tool to use for multiple sequence alignment. Current options are only mafft.

**Default:** `mafft`


### `--download_nextclade_dataset` {#download-nextclade-dataset}

**Type:** `boolean` | *Optional*

Uses included nextclade dataset for SARS-CoV-2 during runtime when false.

**Default:** `True`


### `--predownloaded_nextclade_dataset` {#predownloaded-nextclade-dataset}

**Type:** `string` | *Optional*

Path to predownloaded nextclade dataset


### `--filter` {#filter}

**Type:** `boolean` | *Optional*

Specifies if reference-mapped fastq files should be extracted


### `--markdup` {#markdup}

**Type:** `boolean` | *Optional*

Specifies if duplicate reads should be removed


### `--relatedness` {#relatedness}

**Type:** `boolean` | *Optional*

Turns on multiple sequence alignment subworkflow when true



## Toggles

### `--aci` {#aci}

**Type:** `boolean` | *Optional*

Specifies if process should be used


### `--artic` {#artic}

**Type:** `boolean` | *Optional*

Specifies if process should be used


### `--artic_filter` {#artic-filter}

**Type:** `boolean` | *Optional*

Specifies if process should be used


### `--bbnorm` {#bbnorm}

**Type:** `boolean` | *Optional*

Specifies if process should be used


### `--bcftools_variants` {#bcftools-variants}

**Type:** `boolean` | *Optional*

Specifies if process should be used

**Default:** `True`


### `--fastqc` {#fastqc}

**Type:** `boolean` | *Optional*

Specifies if process should be used

**Default:** `True`


### `--freyja` {#freyja}

**Type:** `boolean` | *Optional*

Specifies if process should be used

**Default:** `True`


### `--freyja_aggregate` {#freyja-aggregate}

**Type:** `boolean` | *Optional*

Specifies if process should be used

**Default:** `True`


### `--heatcluster` {#heatcluster}

**Type:** `boolean` | *Optional*

Specifies if process should be used

**Default:** `True`


### `--igv_reports` {#igv-reports}

**Type:** `boolean` | *Optional*

Specifies if process should be used


### `--iqtree` {#iqtree}

**Type:** `boolean` | *Optional*

Specifies if process should be used

**Default:** `True`


### `--ivar_variants` {#ivar-variants}

**Type:** `boolean` | *Optional*

Specifies if process should be used


### `--kraken2` {#kraken2}

**Type:** `boolean` | *Optional*

Specifies if process should be used

**Default:** `True`


### `--multiqc` {#multiqc}

**Type:** `boolean` | *Optional*

Specifies if process should be used

**Default:** `True`


### `--nextclade` {#nextclade}

**Type:** `boolean` | *Optional*

Specifies if process should be used

**Default:** `True`


### `--pango_aliasor` {#pango-aliasor}

**Type:** `boolean` | *Optional*

Specifies if process should be used

**Default:** `True`


### `--pangolin` {#pangolin}

**Type:** `boolean` | *Optional*

Specifies if process should be used

**Default:** `True`


### `--phytreeviz` {#phytreeviz}

**Type:** `boolean` | *Optional*

Specifies if process should be used

**Default:** `True`


### `--samtools_qc` {#samtools-qc}

**Type:** `boolean` | *Optional*

Specifies if process should be used

**Default:** `True`


### `--samtools_ampliconstats` {#samtools-ampliconstats}

**Type:** `boolean` | *Optional*

Specifies if process should be used

**Default:** `True`


### `--samtools_plot_ampliconstats` {#samtools-plot-ampliconstats}

**Type:** `boolean` | *Optional*

Specifies if process should be used

**Default:** `True`


### `--snpdists` {#snpdists}

**Type:** `boolean` | *Optional*

Specifies if process should be used

**Default:** `True`


### `--vadr` {#vadr}

**Type:** `boolean` | *Optional*

Specifies if process should be used

**Default:** `True`



## Process Adjustments

### `--aci_options` {#aci-options}

**Type:** `string` | *Optional*

Options for process


### `--artic_read_filtering_options` {#artic-read-filtering-options}

**Type:** `string` | *Optional*

Options for process

**Default:** `--min-length 400 --max-length 700`


### `--artic_options` {#artic-options}

**Type:** `string` | *Optional*

Options for process

**Default:** `--normalise 200 --model r1041_e82_400bps_sup_v500 --model-dir /opt/conda/envs/artic/bin/models/`


### `--bbnorm_options` {#bbnorm-options}

**Type:** `string` | *Optional*

Options for process

**Default:** `target=200 min=5`


### `--bcftools_variants_options` {#bcftools-variants-options}

**Type:** `string` | *Optional*

Options for process


### `--fastp_options` {#fastp-options}

**Type:** `string` | *Optional*

Options for process


### `--fastqc_options` {#fastqc-options}

**Type:** `string` | *Optional*

Options for process


### `--filter_options` {#filter-options}

**Type:** `string` | *Optional*

Options for process


### `--heatcluster_options` {#heatcluster-options}

**Type:** `string` | *Optional*

Options for process

**Default:** `-t png`


### `--igv_reports_options` {#igv-reports-options}

**Type:** `string` | *Optional*

Options for process

**Default:** `--flanking 1000`


### `--iqtree_options` {#iqtree-options}

**Type:** `string` | *Optional*

Options for process

**Default:** `-ninit 2 -n 2 -me 0.05 -m GTR`


### `--ivar_consensus_options` {#ivar-consensus-options}

**Type:** `string` | *Optional*

Options for process

**Default:** `-q 20 -t 0.6 -n N`


### `--ivar_trim_options` {#ivar-trim-options}

**Type:** `string` | *Optional*

Options for process


### `--ivar_variants_options` {#ivar-variants-options}

**Type:** `string` | *Optional*

Options for process

**Default:** `-q 20 -t 0.6`


### `--mafft_options` {#mafft-options}

**Type:** `string` | *Optional*

Options for process

**Default:** `--maxambiguous 0.5`


### `--minimap2_options` {#minimap2-options}

**Type:** `string` | *Optional*

Options for process

**Default:** `-K 20M`


### `--multiqc_options` {#multiqc-options}

**Type:** `string` | *Optional*

Options for process


### `--phytreeviz_options` {#phytreeviz-options}

**Type:** `string` | *Optional*

Options for process


### `--pangolin_options` {#pangolin-options}

**Type:** `string` | *Optional*

Options for process


### `--samtools_coverage_options` {#samtools-coverage-options}

**Type:** `string` | *Optional*

Options for process


### `--samtools_ampliconclip_options` {#samtools-ampliconclip-options}

**Type:** `string` | *Optional*

Options for process


### `--samtools_flagstat_options` {#samtools-flagstat-options}

**Type:** `string` | *Optional*

Options for process


### `--samtools_stats_options` {#samtools-stats-options}

**Type:** `string` | *Optional*

Options for process


### `--samtools_ampliconstats_options` {#samtools-ampliconstats-options}

**Type:** `string` | *Optional*

Options for process

**Default:** `--max-amplicon-length 3000 --max-amplicons 3000`


### `--samtools_plot_ampliconstats_options` {#samtools-plot-ampliconstats-options}

**Type:** `string` | *Optional*

Options for process

**Default:** `-size 1200,900 -size2 1200,900 -size3 1200,900`


### `--samtools_depth_options` {#samtools-depth-options}

**Type:** `string` | *Optional*

Options for process


### `--samtools_markdup_options` {#samtools-markdup-options}

**Type:** `string` | *Optional*

Options for process


### `--samtools_fixmate_options` {#samtools-fixmate-options}

**Type:** `string` | *Optional*

Options for process


### `--seqyclean_options` {#seqyclean-options}

**Type:** `string` | *Optional*

Options for process

**Default:** `-minlen 25 -qual`


### `--seqyclean_contaminant_file` {#seqyclean-contaminant-file}

**Type:** `string` | *Optional*

Options for process

**Default:** `/Adapters_plus_PhiX_174.fasta`


### `--snpdists_options` {#snpdists-options}

**Type:** `string` | *Optional*

Options for process

**Default:** `-c`


### `--nextclade_options` {#nextclade-options}

**Type:** `string` | *Optional*

Options for process


### `--vadr_mdir` {#vadr-mdir}

**Type:** `string` | *Optional*

Options for process

**Default:** `/opt/vadr/vadr-models`


### `--freyja_demix_options` {#freyja-demix-options}

**Type:** `string` | *Optional*

Options for process


### `--freyja_update_options` {#freyja-update-options}

**Type:** `string` | *Optional*

Options for process


### `--freyja_variants_options` {#freyja-variants-options}

**Type:** `string` | *Optional*

Options for process


### `--pango_aliasor_options` {#pango-aliasor-options}

**Type:** `string` | *Optional*

Options for process


### `--freyja_aggregate_options` {#freyja-aggregate-options}

**Type:** `string` | *Optional*

Options for process


### `--freyja_plot_options` {#freyja-plot-options}

**Type:** `string` | *Optional*

Options for process


### `--freyja_plot_filetype` {#freyja-plot-filetype}

**Type:** `string` | *Optional*

Options for process

**Default:** `png`


### `--kraken2_options` {#kraken2-options}

**Type:** `string` | *Optional*

Options for process


### `--vadr_options` {#vadr-options}

**Type:** `string` | *Optional*

Options for process


### `--vadr_trim_options` {#vadr-trim-options}

**Type:** `string` | *Optional*

Options for process



## Institutional config options

### `--custom_config_version` {#custom-config-version}

**Type:** `string` | *Optional*

Git commit id for Institutional configs.

**Default:** `master`


### `--custom_config_base` {#custom-config-base}

**Type:** `string` | *Optional*

Base directory for Institutional configs.

> If you're running offline, Nextflow will not be able to fetch the institutional config files from the internet. If you don't need them, then this is not a problem. If you do need them, you should download the files from the repo and tell Nextflow where to find them with this parameter.

**Default:** `https://raw.githubusercontent.com/nf-core/configs/master`


### `--config_profile_name` {#config-profile-name}

**Type:** `string` | *Optional*

Institutional config name.


### `--config_profile_description` {#config-profile-description}

**Type:** `string` | *Optional*

Institutional config description.


### `--config_profile_contact` {#config-profile-contact}

**Type:** `string` | *Optional*

Institutional config contact information.


### `--config_profile_url` {#config-profile-url}

**Type:** `string` | *Optional*

Institutional config URL link.



## Generic options

### `--help` {#help}

**Type:** `boolean` | *Optional*

Display help text.


### `--version` {#version}

**Type:** `boolean` | *Optional*

Display version and exit.


### `--email` {#email}

**Type:** `string` | *Optional*

Email address for completion summary.

> Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

**Pattern:** `^([a-zA-Z0-9_\-\.]+)@([a-zA-Z0-9_\-\.]+)\.([a-zA-Z]{2,5})$`


### `--publish_dir_mode` {#publish-dir-mode}

**Type:** `string` | *Optional*

Method used to save pipeline results to output directory.

> The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.

**Default:** `copy`

**Allowed values:**
- `symlink`
- `rellink`
- `link`
- `copy`
- `copyNoFollow`
- `move`


### `--email_on_fail` {#email-on-fail}

**Type:** `string` | *Optional*

Email address for completion summary, only when pipeline fails.

> An email address to send a summary email to when the pipeline is completed - ONLY sent if the pipeline does not exit successfully.

**Pattern:** `^([a-zA-Z0-9_\-\.]+)@([a-zA-Z0-9_\-\.]+)\.([a-zA-Z]{2,5})$`


### `--plaintext_email` {#plaintext-email}

**Type:** `boolean` | *Optional*

Send plain-text email instead of HTML.


### `--monochrome_logs` {#monochrome-logs}

**Type:** `boolean` | *Optional*

Do not use coloured log outputs.


### `--hook_url` {#hook-url}

**Type:** `string` | *Optional*

Incoming hook URL for messaging service

> Incoming hook URL for messaging service. Currently, MS Teams and Slack are supported.


### `--pipelines_testdata_base_path` {#pipelines-testdata-base-path}

**Type:** `string` | *Optional*

Base URL or local path to location of pipeline test dataset files

**Default:** `https://raw.githubusercontent.com/nf-core/test-datasets/`



---

*This pipeline was built with [Nextflow](https://nextflow.io).
Documentation generated by [nf-docs](https://github.com/ewels/nf-docs) v0.2.0 on 2026-03-12 23:00:43 UTC.*
