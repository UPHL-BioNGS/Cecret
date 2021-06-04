# Cecret

Named after the beautiful [Cecret lake](https://en.wikipedia.org/wiki/Cecret_Lake)

Location: 40.570°N 111.622°W , 9,875 feet (3,010 m) elevation

![alt text](https://upload.wikimedia.org/wikipedia/commons/thumb/c/cb/Cecret_Lake_Panorama_Albion_Basin_Alta_Utah_July_2009.jpg/2560px-Cecret_Lake_Panorama_Albion_Basin_Alta_Utah_July_2009.jpg)

Cecret is a workflow developed by [@erinyoung]("https://github.com/erinyoung") at the [Utah Public Health Laborotory](https://uphl.utah.gov/) for SARS-COV-2 sequencing with the [artic](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html)/Illumina hybrid library prep workflow for MiSeq data with protocols [here](https://www.protocols.io/view/sars-cov-2-sequencing-on-illumina-miseq-using-arti-bffyjjpw) and [here](https://www.protocols.io/view/sars-cov-2-sequencing-on-illumina-miseq-using-arti-bfefjjbn). Built to work on linux-based operating systems. Additional config options are needed for cloud batch usage.



Cecret is also part of the [staphb-toolkit](https://github.com/StaPH-B/staphb_toolkit).

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

### Arrange paired-end fastq.gz reads as follows or designate directory with `params.reads` or `--reads`
```
directory
|-reads
     |-*fastq.gz
```

### Arrange single-end fastq.gz reads as follows or designate directory with `params.single_reads` or `--single_reads`
```
directory
|-single_reads
   |-*fastq.gz
```

WARNING : single and paired-end reads **cannot** be in the same directory

### Start the workflow

```
nextflow run Cecret.nf -c configs/singularity.config
```
## Optional toggles:

### Using fastp to clean reads instead of seqyclean
```
nextflow run Cecret.nf -c configs/singularity.config --cleaner fastp
```
Or set `params.cleaner = 'fastp'` in a config file

### Using samtools to trim amplicons instead of ivar
```
nextflow run Cecret.nf -c configs/singularity.config --trimmer samtools
```
Or set `params.trimmer = 'samtools'` in a config file

### Using minimap2 to align reads instead of bwa
```
nextflow run Cecret.nf -c configs/singularity.config --aligner minimap2
```
Or set `params.aligner = 'minimap2'` in a config file

### Determining relatedness
To create a multiple sequence alignment and corresponding phylogenetic tree and SNP matrix, set `params.relatedness = true` or 
```
nextflow run Cecret.nf -c configs/singularity.config --relatedness true
```
### Classify reads with kraken2
To classify reads with kraken2 to identify reads from human or the organism of your choice
#### Step 1. Get a kraken2 database (note : this link is no longer active. I'm actively working on creating my own)
```
mkdir -p kraken2_db
cd kraken2_db
wget https://storage.googleapis.com/sars-cov-2/kraken2_h%2Bv_20200319.tar.gz
tar -zxf kraken2_h+v_20200319.tar.gz
rm -rf kraken2_h+v_20200319.tar.gz
```
#### Step 2. Set the paramaters accordingly
```
params.kraken2 = true
params.kraken2_db = < path to kraken2 database >
params.kraken2_organism = "Severe acute respiratory syndrome coronavirus 2"
```

## The main components of Cecret are:

- [seqyclean](https://github.com/ibest/seqyclean) - for cleaning reads
- [fastp](https://github.com/OpenGene/fastp) - for cleaning reads ; optional, faster alternative to seqyclean
- [bwa](http://bio-bwa.sourceforge.net/) - for aligning reads to the reference
- [minimap2](https://github.com/lh3/minimap2) - an alternative to bwa
- [ivar](https://andersen-lab.github.io/ivar/html/manualpage.html) - calling variants and creating a consensus fasta; optional primer trimmer
- [samtools](http://www.htslib.org/) - for QC metrics and sorting; optional primer trimmer; optional converting bam to fastq files
- [fastqc](https://github.com/s-andrews/FastQC) - for QC metrics
- [bedtools](https://bedtools.readthedocs.io/en/latest/) - for depth estimation over amplicons
- [kraken2](https://ccb.jhu.edu/software/kraken2/) - for read classification
- [pangolin](https://github.com/cov-lineages/pangolin) - for lineage classification
- [nextclade](https://clades.nextstrain.org/) - for clade classification
- [vadr](https://github.com/ncbi/vadr) - for annotating fastas like NCBI
- [mafft](https://mafft.cbrc.jp/alignment/software/) - for multiple sequence alignment (optional, relatedness must be set to "true")
- [snp-dists](https://github.com/tseemann/snp-dists) - for relatedness determination (optional, relatedness must be set to "true")
- [iqtree](http://www.iqtree.org/) - for phylogenetic tree generation (optional, relatedness must be set to "true")
- [bamsnap](https://github.com/parklab/bamsnap) - to create images of SNPs

### Turning off unneeded processes
It came to my attention that some processes (like bcftools) do not work consistently. Also, they might take longer than wanted and might not even be needed for the end user. Here's the processes that can be turned off with their default values:
```
params.bcftools_variants = false          # the container gets a lot of traffic which can error when attempting to download
params.fastqc = true                      # qc on the sequencing reads
params.ivar_variants = true               # itemize the variants identified by ivar
params.samtools_stats = true              # stats about the bam files
params.samtools_coverage = true           # stats about the bam files
params.samtools_depth = true              # stats about the bam files
params.samtools_flagstat = true           # stats about the bam files
params.samtools_ampliconstats = true      # stats about the amplicons
params.samtools_plot_ampliconstats = true # images related to amplicon performance
params.kraken2 = false                    # used to classify reads and needs a corresponding params.kraken2_db and organism if not SARS-CoV-2
params.bedtools_multicov = true           # bedtools multicov for coverage approximation of amplicons
params.nextclade = true                   # SARS-CoV-2 clade determination
params.pangolin = true                    # SARS-CoV-2 clade determination
params.vadr = false                       # NCBI fasta QC
params.relatedness = false                # actually the toggle for mafft to create a multiple sequence alignement
params.snpdists = true                    # creates snp matrix from mafft multiple sequence alignment
params.iqtree = true                      # creates phylogenetic tree from mafft multiple sequence alignement
params.bamsnap = false                    # can be really slow. Works best with bcftools variants. An example bamsnap image is below.
params.rename = false                     # needs a corresponding sample file and will rename files for GISAID and NCBI submission
params.filter = false                     # takes the aligned reads and turns them back into fastq.gz files
```

### Add Genbank parsable header to consensus fasta

This requires a comma-delimted file set with `params.sample_file` file with a row for each sample and a comma-delimited column for each item to add to the GenBank submission header. Additionally, adjust `params.rename = true`.

The following headers are required
- Sample_id          (required, must match sample_id*.fa*)
- Submission_id      (if file needs renaming)
- Collection_Date

Example covid_samples.csv file contents:
```
Sample_id,Submission_ID,Collection_Date,SRR
12345,UT-UPHL-12345,2020-08-22,SRR1
67890,UT-UPHL-67890,2020-08-22,SRR2
23456,UT-UPHL-23456,2020-08-22,SRR3
78901,UT-UPHL-78901,2020-08-18,SRR4
```
Where the files named `12345-UT-M03999-200822_S9_L001_R1_001.fastq.gz`, `12345-UT-M03999-200822_S9_L001_R2_001.fastq.gz` will be renamed `UT-UPHL-12345.R1.fastq.gz` and `UT-UPHL-12345.R2.fastq.gz`. A GISAID and GenBank friendly multifasta files ready for submission are also generated. The GenBank multifasta uses the input file to create fasta headers like
```
>12345 [Collection_Date=2020-08-22][organism=Severe acute respiratory syndrome coronavirus 2][host=human][country=USA][isolate=SARS-CoV-2/human/USA/12345/2020][SRR=SRR1]
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
```

Sometimes sequencing fails, so there are parameters for how many non-ambiguous bases a fasta needs in order to get incorporated into the final file. This can be set with `params.gisaid_threshold` (Default is '`params.gisaid_threshold = '25000'`') and params.genbank_threshold (Default is '`params.genbank_threshold = '15000'`').

## Final file structure
```
cecret_run_results.txt                # information about the sequencing run that's compatible with legacy workflows
covid_samples.csv                     # only if supplied initially
cecret
|-aligned
| |-pretrimmed.sorted.bam
|-bamsnap
| |-sample
|   |-ivar
|     |-variant.png                   # png of variants identified via ivar
|   |-bcftools
|     |-variant.png                   # png of variants identified via bcftools
|-bedtools_multicov
| |-sample.multicov.txt               # depth per amplicon
|-consensus
| |-sample.consensus.fa               # the likely reason you are running this workflow
|-fastp
| |-sample_clean_PE1.fastq            # clean file: only if params.cleaner=fastp
| |-sample_clean_PE2.fastq            # clean file: only if params.cleaner=fastp
|-fastqc
| |-sample.fastqc.html
| |-sample.fastqc.zip
|-filter                              # optional: turns aligned bams into fastq files
| |-sample_filtered_R1.fastq
| |-sample_filtered_R2.fastq
| |-sample_filtered_unpaired.fastq
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
| |-lineage_report.csv                # identification of pangolin lineages
|-samtools_ampliconstats
| |-sample_ampliconstats.txt          # stats for the amplicons used
|-samtools_coverage
| |-aligned
| | |-sample.cov.hist                 # histogram of coverage for aligned reads
| | |-sample.cov.txt                  # tabular information of coverage for aligned reads
| |-trimmed
|   |-sample.cov.trim.hist            # histogram of coverage for aligned reads after primer trimming
|   |-sample.cov.trim.txt             # tabular information of coverage for aligned reads after primer trimming
|-samtools_depth
| |-aligned
| | |-sample.depth.txt                # read depth for each read position
| |-trimmed
|   |-sample.depth.txt                # read depth for each position
|-samtools_flagstat
| |-aligned
| | |-sample.flagstat.txt             # samtools stats for aligned reads
| |-trimmed
|   |-sample.flagstat.trim.txt        # samtools stats for trimmed reads
|-samtools_plot_ampliconstats
| |-sample.*.png                      # images corresponding to amiplicon performance
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
| |-sample.genbank.fa                 # fasta file with formatting and header including metadata for genbank
| |-sample.gisaid.fa                  # fasta file with header for gisaid
| |-sample.R1.fastq.gz                # renamed raw fastq.gz file
| |-sample.R2.fastq.gz                # renamed raw fastq.gz file
|-summary
| |-sample.summary.txt                # individual results
|-summary.csv                         # tab-delimited summary of results from the workflow
reads
| |-sample_S1_L001_R1_001.fastq.gz    # initial file
| |-sample_S1_L001_R2_001.fastq.gz    # inital file
work                                  # nextflows work directory. Likely fairly large.
vadr
  |-vadr*                             # vadr output
```

# **A FILE THAT THE END USER/YOU CAN COPY AND EDIT IS FOUND AT [configs/example.config](configs/example.config)**
This file contains all of the configurable parameters with their default values. If you are using some sort of cloud or HPC setup, I highly recommend you copy this and edit it appropriately. A limited list of parameters is listed below:

### input and output directories
* params.reads = workflow.launchDir + '/reads'
* params.single_reads = workflow.launchDir + '/single_reads'
* params.outdir = workflow.launchDir + "/cecret"

### reference files for SARS-CoV-2 with artic V3 primers (part of the github repository)
* params.reference_genome = workflow.projectDir + "/configs/MN908947.3.fasta"
* params.gff_file = workflow.projectDir + "/configs/MN908947.3.gff"
* params.primer_bed = workflow.projectDir + "/configs/artic_V3_nCoV-2019.bed"
* params.amplicon_bed = workflow.projectDir + "/configs/nCoV-2019.insert.bed"

### Other useful options
* To "resume" a workflow that use `-resume` with your nextflow command
* To create a report, use `-with-report` with your nextflow command
* To use nextflow tower, use `-with-tower` with your nextflow command

# Frequently Asked Questions (aka FAQ)
### What do I do if I encounter an error?

**TELL ME ABOUT IT!!!**
* [Github issue](https://github.com/UPHL-BioNGS/Cecret/issues)
* Email me
* Send me a message on slack

Be sure to include the command that you used, what config file you used, and what the **nextflow** error was. 

### What if I just want to annotate some SARS-CoV-2 fastas with pangolin, nextclade and vadr?
```
nextflow run Cecret_annotation.nf -c configs/singularity.config --fastas <directory with fastas>
```
You can run mafft, snpdists, and iqtree on a collection of fastas as well with
```
nextflow run Cecret_annotation.nf -c configs/singularity.config --relatedness true --fastas <directory with fastas>
```
[Cecret_annotation.nf](./Cecret_annotation.nf) Uses the same basic structure as the main workflow, so the same config file can probably be used for both.

### Where is an example config file?
You're more than welcome to look at what we use here at UPHL [here](./configs/UPHL.config).

### Why is bcftools set to 'false' by default?

There's nothing wrong with the bcftools process, and the vcf created by bcftools is rather handy for additional analyses. The `'staphb/bcftools:latest'` container is really popular, and has issues downloading during high traffic times. I don't have the time to handle issues of users not understanding why the container did not download. /Sorry

If you want to get the output from bcftools, set `params.bcftools = true` 

### Can I get images of my SNPs and indels?

Yes. Set `params.bamsnap = true`. This is false by default because of how long it takes. It will work with variants called by `ivar` and `bcftools`, although it is **MUCH** faster with the vcf created by bcftools. 

Warning : will not work on all variants. This is due to how bamsnap runs. It is even less likely to work on indels. 

### What is the difference between `params.amplicon_bed` and `params.primer_bed`?

The primer bedfile is the file with the start and stop of each **primer** sequence.

```
$ head -n 3 artic_V3_nCoV-2019.bed 
MN908947.3	30	54	nCoV-2019_1_LEFT	nCoV-2019_1	+
MN908947.3	385	410	nCoV-2019_1_RIGHT	nCoV-2019_1	-
MN908947.3	320	342	nCoV-2019_2_LEFT	nCoV-2019_2	+
```
The amplicon bedfile is the file with the start and stop of each intended **amplicon**.
```
$ head -n 3 nCoV-2019.insert.bed 
MN908947.3	54	385	1	1	+
MN908947.3	342	704	2	2	+
MN908947.3	664	1004	3	1	+
```
Due to the many varieties of primer bedfiles, I determined it was best if the user supplied this file for custom primer sequences.

### What if I am using an amplicon based library that is not SARS-CoV-2?

In your config file, change the following relevant parameters:
```
params.reference_genome
params.primer_bed
params.amplicon_bed or set params.bedtools_multicov = false
params.gff_file or set params.ivar_variants = false
```
And set
```
params.pangolin = false 
params.nextclade = false
params.vadr = false or configure your vadr container appropriately and params.vadr_reference
```
### What if I need to filter out human reads or I only want reads that map to my reference?

Although not perfect, if `params.filter = true`, then only the reads that were mapped to the reference are returned. This should eliminate all human contamination (as long as human is not part of your reference). 

### This workflow has too many bells and whistles. I really only care about generating a consensus fasta. How do I get rid of all the extras?

Change the parameters in your config file and set most of them to false. 

```
params.fastqc = false
params.ivar_variants = false
params.samtools_stats = false
params.samtools_coverage = false
params.samtools_depth = false
params.samtools_flagstat = false
params.bedtools_multicov = false
params.samtools_ampliconstats = false
params.samtools_plot_ampliconstats = false
params.bedtools_multicov = false
params.pangolin = false
params.nextclade = false
params.vadr = false
```
And, yes, this means I added some bells and whistles so you could turn off the bells and whistles. /irony

# Directed Acyclic Diagrams (DAG)
### Full workflow
![alt text](images/Cecret_workflow.png)

# Sample bamsnap plot
![alt text](images/sample_bamsnap.png)

![alt text](https://uphl.utah.gov/wp-content/uploads/New-UPHL-Logo.png)
