# Cecret

Named after the beautiful [Cecret lake](https://en.wikipedia.org/wiki/Cecret_Lake)

Location: 40.570°N 111.622°W , Elevation: 9,875 feet (3,010 m), [Hiking level: easy](https://www.alltrails.com/trail/us/utah/cecret-lake-trail)

|                                                                                                                                                                |                                                                                                                                                                |
| :------------------------------------------------------------------------------------------------------------------------------------------------------------: | :------------------------------------------------------------------------------------------------------------------------------------------------------------: |
| <img src="https://intermountainhealthcare.org/-/media/images/enterpriseserviceslines/live-well/healthy-hikes/cecret-lake/cecret-lake-md-13.jpg" width="500"/>  | <img src="https://intermountainhealthcare.org/-/media/images/enterpriseserviceslines/live-well/healthy-hikes/cecret-lake/cecret-lake-md-1.jpg" width="500" />  |
| <img src="https://intermountainhealthcare.org/-/media/images/enterpriseserviceslines/live-well/healthy-hikes/cecret-lake/cecret-lake-md-10.jpg" width="500" /> | <img src="https://intermountainhealthcare.org/-/media/images/enterpriseserviceslines/live-well/healthy-hikes/cecret-lake/cecret-lake-md-11.jpg" width="500" /> |

([Image credit: Intermountain Healthcare](https://intermountainhealthcare.org/services/wellness-preventive-medicine/live-well/move-well/healthy-hikes/find-a-hike/cecret-lake/))

Table of Contents:

- [Introduction](https://github.com/UPHL-BioNGS/Cecret/wiki/)
- [Dependencies](https://github.com/UPHL-BioNGS/Cecret#dependencies)
- [Usage](https://github.com/UPHL-BioNGS/Cecret/wiki/Usage)
  - [Using a sample sheet](https://github.com/UPHL-BioNGS/Cecret/wiki/Input#using-a-sample-sheet)
- [Input and output directories](https://github.com/UPHL-BioNGS/Cecret/wiki/Input#files-from-directories)
- [Quality Assessment](https://github.com/UPHL-BioNGS/Cecret/wiki/qc)
- [Setting primer and amplicon bedfiles](https://github.com/UPHL-BioNGS/Cecret/wiki/Usage#determining-primer-and-amplicon-bedfiles)
- [Using a predownloaded nextclade dataset](https://github.com/UPHL-BioNGS/Cecret/wiki/Usage#using-the-included-nextclade-dataset)
- [Setting depth for base calls](https://github.com/UPHL-BioNGS/Cecret/wiki/Usage#determining-depth-for-base-calls)
- [SARS-CoV-2 Wastewater](https://github.com/UPHL-BioNGS/Cecret/wiki/sarscov2%E2%80%90wastewater)
- [Monkeypox](https://github.com/UPHL-BioNGS/Cecret/wiki/MPOX)
- [Updating Cecret](https://github.com/UPHL-BioNGS/Cecret#updating-cecret)
- [Optional toggles](https://github.com/UPHL-BioNGS/Cecret/wiki/Usage#optional-toggles)
- [Determining relatedness or creating trees](https://github.com/UPHL-BioNGS/Cecret/wiki/Usage#determining-relatedness)
- [Classified reads with Kraken2](https://github.com/UPHL-BioNGS/Cecret#classify-reads-with-kraken2)
- [Main components](https://github.com/UPHL-BioNGS/Cecret#the-main-components-of-cecret-are)
- [Turning off unneeded processes](https://github.com/UPHL-BioNGS/Cecret#turning-off-unneeded-processes)
- [Final file structure](https://github.com/UPHL-BioNGS/Cecret#final-file-structure)
- [Config files](https://github.com/UPHL-BioNGS/Cecret#config-files)
  - [Using config files](https://github.com/UPHL-BioNGS/Cecret#using-config-files)
- [Frequently Asked Questions (aka FAQ)](https://github.com/UPHL-BioNGS/Cecret#frequently-asked-questions-aka-faq)

## Introduction

Cecret was originally developed by [@erinyoung](https://github.com/erinyoung) at the [Utah Public Health Laborotory](https://uphl.utah.gov/) for SARS-COV-2 sequencing with the [artic](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html)/Illumina hybrid library prep workflow for MiSeq data with protocols [here](https://www.protocols.io/view/sars-cov-2-sequencing-on-illumina-miseq-using-arti-bffyjjpw) and [here](https://www.protocols.io/view/sars-cov-2-sequencing-on-illumina-miseq-using-arti-bfefjjbn). This nextflow workflow, however, is flexible for many additional organisms and primer schemes as long as the reference genome is "_small_" and "_good enough_." In 2022, [@tives82](https://github.com/tives82) added in contributions for Monkeypox virus, including converting IDT's primer scheme to NC_063383.1 coordinates. We are grateful to everyone that has contributed to this repo.

The nextflow workflow was built to work on linux-based operating systems. Additional config options are needed for cloud batch usage.

The library preparation method greatly impacts which bioinformatic tools are recommended for creating a consensus sequence. For example, amplicon-based library prepation methods will required primer trimming and an elevated minimum depth for base-calling. Some bait-derived library prepation methods have a PCR amplification step, and PCR duplicates will need to be removed. This has added complexity and several (admittedly confusing) options to this workflow. Please submit an [issue](https://github.com/UPHL-BioNGS/Cecret/issues) if/when you run into issues.

It is possible to use this workflow to simply annotate fastas generated from any workflow or downloaded from [GISAID](https://www.gisaid.org/) or [NCBI](https://www.ncbi.nlm.nih.gov/sars-cov-2/). There are also options for multiple sequence alignment (MSA) and phylogenetic tree creation from the fasta files.

Cecret is also part of the [staphb-toolkit](https://github.com/StaPH-B/staphb_toolkit).

## Dependencies

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
- [Singularity](https://singularity.lbl.gov/install-linux) or [Docker](https://docs.docker.com/get-docker/) - set the profile as singularity or docker during runtime

## Usage

Cecret can also use a sample sheet for input with the sample name and reads separated by commas. The header must be `sample,fastq_1,fastq_2`. The general rule is the identifier for the file(s), the file locations, and the type if not paired-end fastq files.

Rows match files with their processing needs.

- paired-end reads: `sample,read1.fastq.gz,read2.fastq.gz`
- single-reads reads: `sample,sample.fastq.gz,single`
- nanopore reads : `sample,sample.fastq.gz,ont`
- fasta files: `sample,sample.fasta,fasta`
- multifasta files: `multifasta,multifasta.fasta,multifasta`

Example sample sheet:

```
sample,fastq_1,fastq_2
SRR13957125,/home/eriny/sandbox/test_files/cecret/reads/SRR13957125_1.fastq.gz,/home/eriny/sandbox/test_files/cecret/reads/SRR13957125_2.fastq.gz
SRR13957170,/home/eriny/sandbox/test_files/cecret/reads/SRR13957170_1.fastq.gz,/home/eriny/sandbox/test_files/cecret/reads/SRR13957170_2.fastq.gz
SRR13957177S,/home/eriny/sandbox/test_files/cecret/single_reads/SRR13957177_1.fastq.gz,single
OQ255990.1,/home/eriny/sandbox/test_files/cecret/fastas/OQ255990.1.fasta,fasta
SRR22452244,/home/eriny/sandbox/test_files/cecret/nanopore/SRR22452244.fastq.gz,ont
```

```
# using docker on samples specified in SampleSheet.csv
nextflow run UPHL-BioNGS/Cecret -profile docker --sample_sheet SampleSheet.csv

# using a config file containing all inputs
nextflow run UPHL-BioNGS/Cecret -c file.config
```

Results are roughly organiized into 'params.outdir'/< analysis >/sample.result

A file summarizing all results is found in `'params.outdir'/cecret_results.csv` and `'params.outdir'/cecret_results.txt`.

Consensus sequences can be found in `'params.outdir'/consensus` and end with `*.consensus.fa`.


## Full workflow

![alt text](images/Cecret.png)


## Updating Cecret

```
nextflow pull UPHL-BioNGS/Cecret
```

Cecret has a weekly update schedule. Cecret's versions have three numbers : X.Y.Z. If the first number, X, changes, there has been a major modification. Params may have changed or subworkflows/channels may have been modified. If the second number, Y, changes, there has been a minor to moderate change. These are mainly for bug fixes or the changing the defaults of params. If the last number has been modified, Z, the workflow is basically the same, there have just been some updates in the containers pulled for the workflow. Most of these updates are to keep Freyja, NextClade, and Pangolin current for SARS-CoV-2 analysis.

## The main components of Cecret are:

- [aci](https://github.com/erinyoung/ACI) - for depth estimation over amplicons (optional, set params.aci = true)
- [artic network](https://github.com/artic-network) - for aligning and consensus creation of nanopore reads
- [bbnorm](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbnorm-guide/) - for normalizing reads (optional, set params.bbnorm = true)
- [bcftools](https://samtools.github.io/bcftools/bcftools.html) - for variants
- [bwa](http://bio-bwa.sourceforge.net/) - for aligning reads to the reference
- [fastp](https://github.com/OpenGene/fastp) - for cleaning reads ; (optional, set params.cleaner = 'fastp')
- [fastqc](https://github.com/s-andrews/FastQC) - for QC metrics
- [freyja](https://github.com/andersen-lab/Freyja) - for multiple SARS-CoV-2 lineage classifications
- [heatcluster](https://github.com/DrB-S/heatcluster) - for visualizing SNP matrices generated via SNP dists
- [iqtree](http://www.iqtree.org/) - for phylogenetic tree generation (optional, set params.relatedness = true)
- [igv-reports](https://github.com/igvteam/igv-reports) - visualizing SNPs (optional, set params.igv_reports = true)
- [ivar](https://andersen-lab.github.io/ivar/html/manualpage.html) - calling variants and creating a consensus fasta; default primer trimmer
- [kraken2](https://ccb.jhu.edu/software/kraken2/) - for read classification
- [mafft](https://mafft.cbrc.jp/alignment/software/) - for multiple sequence alignment (optional, relatedness must be set to "true")
- [minimap2](https://github.com/lh3/minimap2) - an alternative to bwa (optional, set params.aligner = minimap2 )
- [multiqc](https://multiqc.info/) - summary of results
- [nextclade](https://clades.nextstrain.org/) - for SARS-CoV-2 clade classification (optional: aligned fasta can be used from this analysis when relatedness is set to "true" and msa is set to "nextclade")
- [pangolin](https://github.com/cov-lineages/pangolin) - for SARS-CoV-2 lineage classification
- [pango aliasor](https://github.com/corneliusroemer/pango_aliasor) - for SARS-CoV-2 lineage tracing
- [phytreeviz](https://github.com/moshi4/phyTreeViz) - for visualizing phylogenetic trees
- [samtools](http://www.htslib.org/) - for QC metrics and sorting; optional primer trimmer; optional converting bam to fastq files; optional duplication marking
- [seqyclean](https://github.com/ibest/seqyclean) - for cleaning reads
- [snp-dists](https://github.com/tseemann/snp-dists) - for relatedness determination (optional, relatedness must be set to "true")
- [vadr](https://github.com/ncbi/vadr) - for annotating fastas like NCBI

### Turning off unneeded processes

It came to my attention that some processes (like bcftools) do not work consistently. Also, they might take longer than wanted and might not even be needed for the end user. Here's the processes that can be turned off with their default values:

```
params.bcftools_variants = true           # vcf of variants
params.fastqc = true                      # qc on the sequencing reads
params.ivar_variants = true               # itemize the variants identified by ivar
params.samtools_stats = true              # stats about the bam files
params.samtools_coverage = true           # stats about the bam files
params.samtools_depth = true              # stats about the bam files
params.samtools_flagstat = true           # stats about the bam files
params.samtools_ampliconstats = true      # stats about the amplicons
params.samtools_plot_ampliconstats = true # images related to amplicon performance
params.kraken2 = false                    # used to classify reads and needs a corresponding params.kraken2_db and organism if not SARS-CoV-2
params.aci = false                        # coverage approximation of amplicons
parms.igv_reports = false                 # SNP IGV images
params.nextclade = true                   # SARS-CoV-2 clade determination
params.pangolin = true                    # SARS-CoV-2 lineage determination
params.pango_aliasor = true               # SARS-CoV-2 lineage tracing
params.freyja = true                      # multiple SARS-CoV-2 lineage determination
params.vadr = false                       # NCBI fasta QC
params.relatedness = false                # create multiple sequence alignments with input fastq and fasta files
params.snpdists = true                    # creates snp matrix from mafft multiple sequence alignment
params.iqtree = true                      # creates phylogenetic tree from mafft multiple sequence alignement
params.bamsnap = false                    # has been removed
params.rename = false                     # needs a corresponding sample file and will rename files for GISAID and NCBI submission
params.filter = false                     # takes the aligned reads and turns them back into fastq.gz files
params.multiqc = true                     # aggregates data into single report
```
