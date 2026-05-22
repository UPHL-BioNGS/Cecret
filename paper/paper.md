# CECRET: A Nextflow pipeline for reference-based amplicon consensus generation

---
title: 'CECRET: A Nextflow pipeline for reference-based amplicon consensus generation'
tags:
  - nextflow
  - bioinformatics
  - genomics
  - public health
  - SARS-CoV-2
  - MPOX
  - amplicon sequencing
authors:
  - name: Erin L. Young
    orcid: 0000-0002-7535-006X
    affiliation: 1
  - name: Thomas Iverson
    orcid: 0009-0002-1670-0928
    affiliation: 1
  - name: John Arnn
    orcid: 0009-0005-2783-3744
    affiliation: 1
  - name: Olinto Linares-Perdomo
    orcid: 0009-0006-4344-7795
    affiliation: 2
  - name: Robert Sainsbury
    orcid: 0009-0001-7253-7777
    affiliation: 1
  - name: Elizabeth Anderson
    orcid: 0000-0002-1050-5470
    affiliation: 1
affiliations:
 - name: Utah Public Health Laboratory, Department of Health and Human Services, State of Utah, Salt Lake City, Utah, USA
   index: 1
 - name: Huntsman Cancer Institute, Department of Oncological Sciences, University of Utah, Salt Lake City, Utah, USA
   index: 2
date: 24 April 2026
bibliography: paper.bib
---

## Summary

CECRET is an open-source bioinformatics workflow for reference-based consensus genome generation from amplicon sequencing data. Developed at the Utah Public Health Laboratory (UPHL), the workflow is optimized for Illumina sequencing using ARTIC-style primer schemes and was designed to support high-throughput viral genomic surveillance in public health laboratories. The pipeline automates preprocessing, alignment, primer trimming, consensus generation, quality assessment, and lineage assignment to produce standardized analysis-ready outputs from raw sequencing reads.

Unlike many early viral surveillance workflows that were primarily designed for Oxford Nanopore Technologies (ONT) data,[@Quick:2017] CECRET emphasizes accurate primer trimming, depth-aware filtering, and reproducible consensus generation for amplicon-based Illumina sequencing. The workflow is implemented in Nextflow with containerized dependencies,[@DiTommaso:2017; @Kurtzer:2017; @Florek:2025] enabling reproducible execution across local, high-performance computing, cloud, and restricted-access environments. Although initially developed for SARS-CoV-2 surveillance, CECRET has also been adapted for pathogens including MPOX and Measles virus. The software is released under the MIT license and is publicly available through GitHub and Zenodo [@cecret].


## Statement of Need

Amplicon-based sequencing protocols derived from the ARTIC framework became widely adopted during the SARS-CoV-2 pandemic because they enabled rapid and cost-effective viral genomic surveillance [@Quick:2017]. Many public health laboratories implemented these protocols using Illumina sequencing platforms [@Bull:2020; @Grubaugh:2019]. However, early bioinformatics workflows were primarily optimized for Oxford Nanopore Technologies (ONT) data and did not consistently enforce processing steps required for reliable Illumina amplicon consensus generation, including post-alignment primer trimming and depth-aware consensus filtering.

CECRET was developed to provide a standardized, reference-based workflow specifically tailored to amplicon-based Illumina sequencing in public health surveillance settings. The workflow encodes operational best practices for preprocessing, alignment, primer trimming, consensus generation, and lineage assignment while minimizing manual parameter tuning and local customization. This design emphasizes reproducibility and consistency across laboratories generating large volumes of surveillance data.

General-purpose workflows such as nf-core/viralrecon [@viralrecon:2024] support diverse sequencing strategies including metagenomics and de novo assembly. In contrast, CECRET prioritizes a constrained, production-oriented design optimized for high-throughput viral surveillance using small reference-guided genomes. The workflow is intended for public health laboratories, clinical genomics groups, and researchers requiring rapid and standardized consensus genome generation from amplicon sequencing data.

## State of the Field

Several workflows are available for viral genome reconstruction and genomic surveillance, including the ARTIC bioinformatics pipelines [@Quick:2017] and more general frameworks such as nf-core/viralrecon [@viralrecon:2024]. These workflows provide broad support for diverse sequencing platforms, metagenomic analyses, and flexible research-oriented configurations.

CECRET was developed to address a narrower operational niche: high-throughput consensus genome generation from amplicon-based Illumina sequencing in public health laboratory settings. In contrast to general-purpose workflows that prioritize flexibility across sequencing strategies, CECRET emphasizes constrained, reference-guided analysis with standardized defaults and minimal configuration requirements. The workflow encodes practices commonly required for Illumina amplicon surveillance, including post-alignment primer trimming, depth-aware consensus generation, and integrated lineage assignment.

Developing CECRET as a standalone workflow rather than contributing to an existing general-purpose framework allowed these operational assumptions and public health reporting requirements to be implemented as core workflow behaviors rather than optional configuration layers. This design prioritizes reproducibility, portability, and consistency across laboratories processing large volumes of surveillance samples under routine production conditions.

## Build vs. Contribute Justification

CECRET was developed as a standalone workflow because existing pipelines prioritized either Oxford Nanopore Technologies (ONT) sequencing or broad multi-platform flexibility rather than operational standardization for Illumina amplicon surveillance. Public health laboratories required a workflow that enforced domain-specific best practices — including post-alignment primer trimming, depth-aware consensus generation, and standardized lineage assignment — as default behaviors rather than optional configuration choices.

Implementing these assumptions within a dedicated workflow enabled a constrained, production-oriented design optimized for reproducibility, minimal parameter tuning, and consistent outputs across laboratories. In contrast, general-purpose frameworks such as nf-core/viralrecon [@viralrecon:2024] are designed to support a broader range of sequencing strategies and analytical objectives.

This design philosophy also enabled tighter integration of pathogen-specific surveillance utilities, quality-control metrics, and operation in restricted or air-gapped public health computing environments.

## Software Design

CECRET is implemented in Nextflow DSL2 [@DiTommaso:2017], enabling modular workflow composition, reproducible execution, and portability across local systems, high-performance computing clusters, and cloud environments. The workflow is organized into subworkflows representing major stages of amplicon-based consensus generation, including preprocessing, alignment, primer trimming, consensus generation, quality assessment, and downstream lineage analysis.

A central design goal of CECRET was operational standardization for public health surveillance workflows. Rather than maximizing configurability across sequencing strategies, the pipeline adopts a constrained reference-guided architecture optimized for Illumina amplicon sequencing of small viral genomes. This design reduces the need for local customization and enforces consistent analytical behaviors across laboratories and sequencing runs.

The workflow integrates established bioinformatics tools for alignment, primer trimming, consensus generation, and lineage assignment, while exposing only a limited set of operationally relevant parameters. Depth-aware filtering and post-alignment primer trimming are implemented as core workflow behaviors to reduce amplification artifacts and improve reproducibility of consensus genome generation from tiled amplicon sequencing data.

All dependencies are distributed through containerized environments using Docker, Singularity, or Apptainer [@Kurtzer:2017], minimizing variability caused by local system configuration. Predefined execution profiles support local, Slurm, and containerized execution environments, while continuous integration testing and parameter validation improve reproducibility and workflow stability.

Because public health laboratories frequently operate in restricted-access or air-gapped computing environments, CECRET was designed to minimize runtime dependence on external network resources. The workflow supports predownloaded reference datasets, bundled lineage assignment resources, and optional disabling of database updates to ensure stable execution in secure environments.

Optional downstream modules support phylogenetic analysis, SNP distance calculation, and wastewater lineage deconvolution, allowing the workflow to extend beyond consensus generation while preserving the same reproducible execution framework.

## Research Impact Statement

CECRET has been used operationally for viral genomic surveillance in public health laboratories during the SARS-CoV-2 pandemic and subsequent outbreak responses. At the Utah Public Health Laboratory (UPHL), the workflow supported routine generation of consensus genomes and lineage assignments from Illumina amplicon sequencing data contributing to statewide surveillance efforts (BioProject PRJNA614995). The workflow has subsequently been adapted for additional pathogens including MPOX (PRJNA843095) and Measles virus (PRJNA1293457).

The software is distributed through the StaPH-B toolkit [@Florek:staphb], supporting adoption across public health laboratory networks and integration within a broader ecosystem of reproducible bioinformatics workflows. The repository includes continuous integration testing, predefined execution profiles, and reproducible test datasets to support validation across computing environments.

Within the StaPH-B consortium, CECRET is widely used for SARS-CoV-2 consensus genome generation from Illumina amplicon sequencing data [@CDCgov:SARS2Sequencing]. The CDC Enteric Diseases Laboratory Branch, with support from the Respiratory Viruses Branch, developed SC2CLIA Cecret, a CLIA-ready derivative workflow that extends CECRET with CDC-specific quality-control metrics, reporting, and automated NCBI submission capabilities [@CDCgov:SC2CLIA].

Together, these deployments demonstrate the workflow’s utility as a reproducible and production-oriented platform for high-throughput viral genomic surveillance in operational public health settings.

## Validation and Benchmarking

CECRET was evaluated on 21 amplicon sequencing samples (15 SARS-CoV-2 and 6 Measles virus) and compared against nf-core/viralrecon. Consensus genomes generated by both workflows were aligned to common reference genomes and compared using pairwise SNP differences and lineage assignment concordance.

Across the 19 samples successfully processed by both workflows, consensus sequences differed by a median of 0 SNPs (range: 0–7). CECRET successfully generated consensus genomes for all 21 samples, while nf-core/viralrecon failed on two low-coverage samples (<100,000 reads). Lineage assignment concordance was 76.9%, with discrepancies primarily occurring in low-quality or high-N consensus genomes classified as “Unassigned.” Measles clade assignment showed 100% concordance across all samples.

Overall, these results indicate that CECRET produces consensus genomes consistent with an established viral genomics workflow while maintaining robust performance in low-coverage conditions commonly encountered in routine surveillance datasets.

## AI Usage Disclosure

Generative AI tools (ChatGPT and Google Gemini) were used to assist with debugging, workflow refactoring, documentation drafting, and manuscript outlining. All generated code, text, and architectural suggestions were manually reviewed and validated by the maintainer for technical accuracy and biological relevance.
