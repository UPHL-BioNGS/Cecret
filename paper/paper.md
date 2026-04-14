---
title: 'Cecret: A Nextflow pipeline for reference-based amplicon consensus generation'
tags:
  - nextflow
  - bioinformatics
  - genomics
  - public health
  - SARS-CoV-2
  - Mpox
  - amplicon sequencing
authors:
  - name: Erin Young
    orcid: 0000-0002-7535-006X
    affiliation: 1
    # Add other authors here
affiliations:
 - name: Utah Public Health Laboratory, Department of Health and Human Services, State of Utah, Salt Lake City, Utah, USA
   index: 1
date: 14 April 2026
bibliography: paper.bib
---

# Summary

Cecret is an open-source Nextflow workflow designed for reference-based amplicon consensus generation. Originally developed for SARS-CoV-2 sequencing using the ARTIC network protocols and Illumina/Nanopore hybrid library workflows, Cecret has been expanded to support a variety of organisms and primer schemes, including Mpox (Monkeypox virus) and other viral pathogens. The pipeline streamlines the complex bioinformatic processes required to go from raw sequencing reads to high-quality consensus sequences, lineage classifications, and phylogenetic trees.

# Statement of need

Public health laboratories and genomic surveillance teams frequently utilize amplicon-based sequencing to monitor viral outbreaks. However, generating accurate consensus sequences from these protocols requires managing complex, multi-step bioinformatics workflows that include adapter trimming, primer clipping, variant calling, and lineage assignment. The specific library preparation method chosen greatly impacts which bioinformatic tools are required.

Cecret addresses this need by providing a flexible, reproducible, and highly configurable Nextflow pipeline [@DiTommaso:2017]. It integrates established bioinformatics tools into a containerized workflow (Docker/Singularity/Apptainer) that can be executed on local machines, HPC clusters, or cloud infrastructure. Cecret handles both paired-end (Illumina) and single-end (Nanopore) reads and incorporates organism-specific typing tools such as Pangolin [@OToole:2021], Nextclade [@Aksamentov:2021], and Freyja [@Karthikeyan:2022]. Furthermore, Cecret handles downstream phylogenetic analysis using MAFFT and IQ-TREE, and aggregates quality control metrics via MultiQC, making it an essential, production-ready tool for genomic epidemiology.

# Features and Functionality

The pipeline accepts standard FASTQ files (or downloaded accessions) and performs:
1. **Read Cleaning and QC:** utilizing Seqyclean or fastp, followed by FastQC.
2. **Alignment:** mapping reads to a reference genome using BWA or minimap2.
3. **Primer Trimming:** securely clipping amplicon primers using iVar or Samtools.
4. **Consensus Generation and Variant Calling:** utilizing iVar, BCFtools, or the ARTIC pipeline.
5. **Taxonomic Classification:** using Kraken2 to detect contamination.
6. **Lineage and Clade Assignment:** integrating Freyja, Pangolin, and Nextclade for rapid pathogen typing.

# Acknowledgements

We acknowledge contributions from the Utah Public Health Laboratory and community contributors who expanded the pipeline's capabilities for Mpox and other pathogens. 

# References
