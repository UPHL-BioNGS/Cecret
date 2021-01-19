# Cecret

Named after the beautiful [Cecret lake](https://en.wikipedia.org/wiki/Cecret_Lake)

Location: 40.570°N 111.622°W , 9,875 feet (3,010 m) elevation

![alt text](https://upload.wikimedia.org/wikipedia/commons/thumb/c/cb/Cecret_Lake_Panorama_Albion_Basin_Alta_Utah_July_2009.jpg/2560px-Cecret_Lake_Panorama_Albion_Basin_Alta_Utah_July_2009.jpg)

Cecret is a workflow developed by @erinyoung at the [Utah Public Health Laborotory](https://uphl.utah.gov/) for SARS-COV-2 sequencing with the [artic](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html)/Illumina hybrid library prep workflow for MiSeq data with protocols [here](https://www.protocols.io/view/sars-cov-2-sequencing-on-illumina-miseq-using-arti-bffyjjpw) and [here](https://www.protocols.io/view/sars-cov-2-sequencing-on-illumina-miseq-using-arti-bfefjjbn). Built to work on linux-based operating systems. Additional config files are needed for cloud batch usage.

# Getting started

```
git clone https://github.com/UPHL-BioNGS/Cecret.git
```

To make your life easier, follow with

```
cd Cecret
git init
```

so that you can use `git pull` for updates.

## Prior to starting the workflow

### Install dependencies
- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
   - Nextflow version 20+ is required (`nextflow -v` to check your installation)
- [Singularity](https://singularity.lbl.gov/install-linux)

or
- [Docker](https://docs.docker.com/get-docker/) (*with the caveat that the creator and maintainer uses singularity and may not be able to troubleshoot all docker issues*)

# Usage

### Arrange paired-end fastq.gz reads as follows or designate directory with `--reads`
```
directory
|-Sequencing_reads
   |-Raw
     |-*fastq.gz
```

### Arrange single-end fastq.gz reads as follows or designate directory with `--single_reads`
```
directory
|-Sequencing_reads
   |-Single
     |-*fastq.gz
```

WARNING : single and paired-end reads **cannot** be in the same directory

### Start the workflow

```
nextflow run Cecret.nf -c config/singularity.config
```

## Optional usage:

### Determining relatedness
Resuming the workflow to create a multiple sequence alignment, phylogenetic tree, and count SNPs by setting the relatedness paramter to `true`
```
nextflow run Cecret.nf -c config/singularity.config -resume --relatedness true
```

### Using samtools to trim amplicons instead of ivar
Setting trimmer parameter to `samtools`
```
nextflow run Cecret.nf -c config/singularity.config -resume --trimmer samtools
```

### Using fastp to clean reads instead of seqyclean
Setting cleaner parameter to `fastp`
```
nextflow run Cecret.nf -c config/singularity.config -resume --cleaner fastp
```

### Classify reads with kraken2
Example download of kraken2 database for human and viral reads (including SARS-CoV-2)
```
mkdir -p kraken2_db
cd kraken2_db
wget https://storage.googleapis.com/sars-cov-2/kraken2_h%2Bv_20200319.tar.gz
tar -zxf kraken2_h+v_20200319.tar.gz
rm -rf kraken2_h+v_20200319.tar.gz
```

```
nextflow run Cecret.nf -c config/singularity.config --kraken2 true --kraken2_db=kraken2_db
```

### Add Genbank parsable header to consensus fasta

This requires a covid_samples.txt file with a row for each sample and a tab-delimited column for each item to add to the header.

The following headers are accepted
- Sample_id          (required, must match sample_id*.fa*)
- Submission_id      (if file needs renaming)
- Country            (default is 'USA')
- Host               (default is 'Human')
- Isolate            (default is submission_id)
- Collection_Date
- Isolation_Source   (default is 'SARS-CoV-2/host/location/isolate/date')
- Clone
- Collected_By
- Fwd_Primer_Name
- Fwd_Primer_Seq
- Latitude_Longitude
- Rev_Primer_Name
- Rev_Primer_Seq
- Note
- Bioproject
- Biosample
- Sra

Example covid_samples.txt file contents:
```
Lab_Accession	Submission_ID	Collection_Date	SRR
12345	UT-UPHL-12345	2020-08-22	SRR1
67890	UT-UPHL-67890	2020-08-22	SRR2
23456	UT-UPHL-23456	2020-08-22	SRR3
78901	UT-UPHL-78901	2020-08-18	SRR4
```
Where the files named `12345-UT-M03999-200822_S9_L001_R1_001.fastq.gz`, `12345-UT-M03999-200822_S9_L001_R2_001.fastq.gz` will be renamed `UT-UPHL-12345.R1.fastq.gz` and `UT-UPHL-12345.R2.fastq.gz`. A consensus file will be duplicated and named `UT-UPHL-12345.consensus.fa`. GISAID and GenBank friendly multifasta files ready for submission are also generated. The GenBank multifasta uses the input file to create fasta headers like
```
>12345 [Collection_Date=2020-08-22][organism=Severe acute respiratory syndrome coronavirus 2][host=human][country=USA][isolate=SARS-CoV-2/human/USA/12345/2020][SRR=SRR1]
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
```
### Turning off unneeded processes
It came to my attention that some processes (like bcftools) do not work consistently. Also, they might take longer than wanted and might not even be needed for the end user. Here's the processes that can be turned off with their default values:
```
params.prepare_reference = true
params.bcftools_variants = false      # the container gets a lot of traffic which can error when attempting to download
params.fastqc = true
params.ivar_variants = true
params.samtools_stats = true
params.samtools_coverage = true
params.samtools_flagstat = true
params.samtools_ampliconstats = true
params.kraken2 = false                # needs a corresponding params.kraken2_db
params.bedtools = true
params.nextclade = true
params.pangolin = true
params.bamsnap = false                # doesn't actually work
params.relatedness = false            # actually the toggle for mafft
params.snpdists = true
params.iqtree = true
```

## The main components of Cecret are:

- [seqyclean](https://github.com/ibest/seqyclean) - for cleaning reads
- [fastp](https://github.com/OpenGene/fastp) - for cleaning reads ; optional, faster alternative to seqyclean
- [bwa](http://bio-bwa.sourceforge.net/) - for aligning reads to the reference
- [minimap2](https://github.com/lh3/minimap2) - an alternative to bwa
- [ivar](https://andersen-lab.github.io/ivar/html/manualpage.html) - calling variants and creating a consensus fasta; optional primer trimmer
- [samtools](http://www.htslib.org/) - for QC metrics and sorting; optional primer trimmer
- [fastqc](https://github.com/s-andrews/FastQC) - for QC metrics
- [bedtools](https://bedtools.readthedocs.io/en/latest/) - for depth estimation over amplicons
- [kraken2](https://ccb.jhu.edu/software/kraken2/) - for read classification
- [pangolin](https://github.com/cov-lineages/pangolin) - for lineage classification
- [nextclade](https://clades.nextstrain.org/) - for clade classification
- [mafft](https://mafft.cbrc.jp/alignment/software/) - for multiple sequence alignment (optional, relatedness must be set to "true")
- [snp-dists](https://github.com/tseemann/snp-dists) - for relatedness determination (optional, relatedness must be set to "true")
- [iqtree](http://www.iqtree.org/) - for phylogenetic tree generation (optional, relatedness must be set to "true")


## Final file structure
```
run_results.txt                       # information about the sequencing run that's compatible with legacy workflows
covid_samples.txt                     # only if supplied initially
cecret
|-aligned
| |-pretrimmed.sorted.bam
|-bedtools
| |-sample.multicov.txt               # depth per amplicon
|-consensus
| |-consensus.fa                      # the likely reason you are running this workflow
|-fastp
| |-sample_clean_PE1.fastq            # clean file: only if params.cleaner=fastp
| |-sample_clean_PE2.fastq            # clean file: only if params.cleaner=fastp
|-fastqc
| |-sample.fastqc.html
| |-sample.fastqc.zip
|-iqtree                              # optional: relatedness parameter must be set to true
| |-iqtree.treefile
|-ivar_trim
| |-sample.primertrim.bam             # aligned reads after primer trimming. trimmer parameter must be set to 'ivar'
|-ivar_variants
| |-sample.variants.tsv               # list of variants identified via ivar and corresponding scores
|-kraken2
| |-sample_kraken2_report.txt         # kraken2 report of the percentage of reads matching virus and human sequences
|-logs
| |-process_logs                      # for troubleshooting puroses
|-mafft                               # optional: relatedness parameter must be set to true
| |-mafft_aligned.fasta               # multiple sequence alignement generated via mafft
|-nextclade                           # identfication of nextclade clades and variants identified
| |-sample_nextclade_report.csv       # actually a ";" deliminated file
|-pangolin
| |-sample
|   |-lineage_report.csv              # identification of pangolin lineages
|-samtools_coverage
| |-aligned
| | |-sample.cov.hist                 # histogram of coverage for aligned reads
| | |-sample.cov.txt                  # tabular information of coverage for aligned reads
| |-trimmed
|   |-sample.cov.trim.hist            # histogram of coverage for aligned reads after primer trimming
|   |-sample.cov.trim.txt             # tabular information of coverage for aligned reads after primer trimming
|-samtools_flagstat
| |-aligned
| | |-sample.flagstat.txt             # samtools stats for aligned reads
| |-trimmed
|   |-sample.flagstat.trim.txt        # samtools stats for trimmed reads
|-samtools_stats
| |-aligned
| | |-sample.stats.txt                # samtools stats for aligned reads
| |-trimmed
|   |-sample.stats.trim.txt           # samtools stats for trimmed reads
|-seqyclean
| |-sample_clean_PE1.fastq            # clean file
| |-sample_clean_PE2.fastq            # clean file
|-snp-dists                           # optional: relatedness parameter must be set to true
| |-snp-dists                         # file containing a table of the number of snps that differ between any two samples
|-submission_files                    # optional: is only created if covid_samples.txt exists
| |-cecret.genbank_submission.fasta   # multifasta for direct genbank
| |-cecret.gisaid_submission.fasta    # multifasta file bulk upload to gisaid
| |-sample.consensus.fa               # renamed consensus fasta file
| |-sample.genbank.fa                 # fasta file with formatting and header including metadata for genbank
| |-sample.gisaid.fa                  # fasta file with header for gisaid
| |-sample.R1.fastq.gz                # renamed raw fastq.gz file
| |-sample.R2.fastq.gz                # renamed raw fastq.gz file
|-summary
| |-sample.summary.txt                # individual results
|-summary.txt                         # tab-delimited summary of results from the workflow
Sequencing_reads
|-Raw
| |-sample_S1_L001_R1_001.fastq.gz    # initial file
| |-sample_S1_L001_R2_001.fastq.gz    # inital file
work                                  # nextflows work directory. Likely fairly large.
```

# Adjustable Paramters with their default values
Parameters can be adjusted in a config file or on the command line. Command line adjustments look like --trimmer 'samtools'

### input and output directories
* params.reads = workflow.launchDir + '/Sequencing_reads/Raw'
* params.single_reads = workflow.launchDir + '/Sequencing_reads/Single'
* params.outdir = workflow.launchDir + "/cecret"
* params.sample_file = workflow.launchDir + '/covid_samples.txt' (optional)

### reference files for SARS-CoV-2 (part of the github repository)
* params.reference_genome = workflow.projectDir + "/config/MN908947.3.fasta"
* params.gff_file = workflow.projectDir + "/config/MN908947.3.gff"
* params.primer_bed = workflow.projectDir + "/config/artic_V3_nCoV-2019.bed"

### toggles for trimmers and cleaner
* params.trimmer = 'ivar'
* params.cleaner = 'seqyclean'
* params.aligner = 'bwa'

### for minimap2
* params.minimap2_K = '20M'

### for ivar trimming and variant/consensus calling
* params.ivar_quality = 20
* params.ivar_frequencing_threshold = 0.6
* params.ivar_minimum_read_depth = 10
* params.mpileup_depth = 8000

### for optional kraken2 contamination detection
* params.kraken2 = false
* params.kraken2_db = ''

### for optional route of tree generation and counting snps between samples
* params.relatedness = false
* params.max_ambiguous = '0.50'
* params.outgroup = 'MN908947.3'
* params.mode='GTR'

### CPUS to use
* params.maxcpus = Runtime.runtime.availableProcessors()

### for file submission headers
* params.year = Calendar.getInstance().get(Calendar.YEAR)
* params.country = 'USA'

### contaminant file for seqyclean (inside the seqyclean container)
* params.seqyclean_contaminant_file="/Adapters_plus_PhiX_174.fasta"
* params.seqyclean_minlen = 25

### Other useful options
To create a report, use `-with-report` with your nextflow command.

# Directed Acyclic Diagrams (DAG)
### Full workflow
![alt text](images/Cecret_workflow.png)

![alt text](https://uphl.utah.gov/wp-content/uploads/New-UPHL-Logo.png)
