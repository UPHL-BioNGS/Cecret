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
  - name: Olinto Linares-Perdomo
    orcid: 0009-0006-4344-7795
    affiliation: 2
  - name: John Arnn
    orcid: 0009-0005-2783-3744
    affiliation: 1
    # Add other authors here
affiliations:
 - name: Utah Public Health Laboratory, Department of Health and Human Services, State of Utah, Salt Lake City, Utah, USA
   index: 1
 - name: Huntsman Cancer Institute, Department of Oncological Sciences, University of Utah, Salt Lake City, Utah, USA
   index: 2
date: 24 April 2026
bibliography: paper.bib
---

## Summary

CECRET is an open-source bioinformatics workflow for reference-based consensus genome generation from amplicon sequencing data. The pipeline was developed at the Utah Public Health Laboratory (UPHL) to address a specific operational need of producing high-quality consensus sequences from Illumina data generated using ARTIC-style amplicon protocols. CECRET is released under the MIT license The source code and archived release are available at https://github.com/UPHL-BioNGS/Cecret and Zenodo (DOI: 10.5281/zenodo.16414959) [@cecret].

While widely adopted bioinformatics workflows exist, many were initially optimized for Oxford Nanopore Technologies (ONT) data [@Quick:2017]. In contrast, CECRET is purpose-built for amplicon-based Illumina sequencing, where accurate primer trimming, depth-aware filtering, and standardized downstream lineage assignment are critical for reliable consensus generation.

The workflow automates the transformation of raw sequencing reads into analysis-ready outputs through a modular pipeline that includes read preprocessing, reference alignment, primer trimming, consensus generation, and downstream pathogen typing. It integrates established tools such as BWA or minimap2 for alignment [@Li:2009; @Li:2018], iVar for primer trimming and consensus generation [@Grubaugh:2019], and lineage and clade assignment tools including Nextclade and Pangolin [@Aksamentov:2021; @OToole:2021]. Optional modules support multiple sequence alignment and phylogenetic inference for outbreak investigation.

CECRET is implemented in Nextflow and distributed with containerized dependencies (Docker, Singularity, or Apptainer) utilizing the State Public Health Bioinformatics community’s containerized software repository [@DiTommaso:2017; @Kurtzer:2017; @Florek:2025], enabling reproducible execution across local, high-performance computing, and cloud environments. Although initially developed for SARS-CoV-2, the workflow has been extended to support additional viral pathogens such as MPOX and Measles Virus, and can be adapted to other small viral genomes given appropriate reference and primer scheme inputs.

## Statement of Need

Rapid and reproducible genomic surveillance is essential for monitoring viral evolution and informing public health responses during outbreaks [@Tiwari:2025; @Florek:2025]. Amplicon-based sequencing approaches, particularly those derived from ARTIC protocols, have been widely adopted due to their scalability and cost efficiency. However, the bioinformatic processing of these data presents distinct challenges that are not consistently addressed by general-purpose workflows.

At the onset of the SARS-CoV-2 pandemic, existing bioinformatics pipelines provided by the ARTIC Network were primarily optimized for Oxford Nanopore Technologies (ONT) data [@Quick:2017]. Public health laboratories frequently implemented hybrid workflows combining ARTIC primer schemes with Illumina sequencing platforms [@Bull:2020; @Grubaugh:2019], creating a gap in available tools. In these settings, accurate consensus generation requires post-alignment primer trimming and depth-aware variant filtering to avoid technical artifacts introduced during amplification.

CECRET was developed to address this gap by providing a standardized, reference-based workflow tailored specifically to amplicon-based Illumina data. The pipeline encodes best practices for:

* post-alignment primer trimming to remove amplification artifacts,
* minimum depth thresholds for high-confidence consensus generation,
* and automated integration of lineage and clade assignment tools.

In contrast to general-purpose workflows such as nf-core/viralrecon (e.g., nf-core/viralrecon [@viralrecon:2024]), which support diverse sequencing strategies including shotgun metagenomics and *de novo* assembly, CECRET prioritizes a constrained, reference-guided design optimized for high-throughput public health surveillance. This design enables consistent, reproducible outputs with minimal parameter tuning, making it suitable for routine operational use.

The primary users of CECRET are public health laboratories, clinical genomics groups, and researchers conducting targeted viral surveillance using amplicon sequencing. The workflow is particularly suited to scenarios requiring rapid turnaround and standardized outputs, such as outbreak response and longitudinal monitoring programs.

## State of the Field

Early in the SARS-CoV-2 pandemic, viral genomic surveillance workflows were often assembled from heterogeneous tools with limited standardization. While the ARTIC Network provided widely adopted laboratory protocols and associated bioinformatics guidance, these pipelines were primarily designed for ONT data and required adaptation for other sequencing platforms [@Quick:2017]. Public health laboratories using Illumina-based amplicon sequencing frequently relied on custom scripts and locally maintained workflows, resulting in variability in analytical results and reduced reproducibility across institutions.

Since 2020, the field has progressed toward standardized, containerized workflows that emphasize reproducibility and portability. Frameworks such as Nextflow and nf-core have enabled the development of community-maintained pipelines that support a wide range of sequencing modalities and analytical goals [@Langer:2025]. These include comprehensive workflows capable of handling metagenomic data, reference-based analyses, and *de novo* assembly within a single framework.

Despite these advances, there remains a need for workflows that prioritize operational simplicity and encode domain-specific best practices for targeted surveillance applications. In particular, amplicon-based Illumina sequencing introduces methodological constraints, such as the requirement for post-alignment primer trimming, that are not consistently enforced in general-purpose pipelines.

CECRET addresses this niche by providing a focused, reference-based workflow optimized for small viral genomes and high-throughput surveillance contexts. Rather than maximizing flexibility across sequencing strategies, the pipeline emphasizes robustness, reproducibility, and minimal configuration for a well-defined class of analyses commonly performed in public health laboratories.

## Build vs. Contribute Justification

CECRET was developed as a standalone workflow rather than as a contribution to an existing pipeline due to the specific operational constraints of public health genomic surveillance during the early stages of the SARS-CoV-2 pandemic.

At that time, available workflows were either narrowly tailored to ONT data (e.g., ARTIC pipelines) or designed as flexible research frameworks supporting multiple sequencing strategies. These approaches did not provide a stable, Illumina-focused solution with enforced best practices for amplicon data processing. In particular, critical steps such as post-alignment primer trimming and depth-based consensus filtering were not consistently implemented as default behaviors.

Developing a dedicated workflow enabled the encoding of these domain-specific requirements as fixed or strongly guided defaults, reducing the need for local customization and minimizing variability between laboratories. This design approach prioritizes reproducibility and consistency over configurability, which is a key requirement in public health settings where standardized outputs are necessary for downstream reporting and data sharing.

In addition, CECRET integrates pathogen-specific downstream tools (e.g., lineage and clade assignment, wastewater deconvolution) into a single pipeline, reducing the need for manual chaining of independent tools. This level of integration would have required substantial restructuring of existing general-purpose workflows.

As the field has matured, comprehensive pipelines such as nf-core/viralrecon has expanded in scope and capability. However, CECRET continues to serve a complementary role by providing a constrained, production-oriented workflow optimized for high-throughput amplicon-based surveillance using Illumina data.

## Software Design

CECRET is implemented using the Nextflow workflow manager (DSL2), enabling modular composition, reproducibility, and portability across diverse computational environments [@DiTommaso:2017]. The pipeline is organized into discrete subworkflows and processes that reflect the logical stages of amplicon-based consensus generation.

### Architecture and Modularity

The workflow is structured into several primary subworkflows:

* Initialization: Standardizes input formats (paired-end Illumina reads, single-end reads, nanopore reads, and FASTA inputs), validates parameters, and prepares reference genome resources.
* Consensus Generation: Performs read preprocessing, alignment, and primer trimming. Preprocessing consists of filtering Illumina reads with either SeqyClean or fastp, followed by an optional normalization step using BBNorm [@Zhbannikov:2017; @Chen:2018; @Bushnell:2014]. Illumina reads are then aligned using BWA by default (with optional minimap2 support), followed by primer trimming and consensus generation using iVar [@Grubaugh:2019], although samtools ampliconclip is also supported. For nanopore data, ARTIC-based tools are used for filtering and consensus generation.
* Quality Control and Variant Calling: Generates coverage metrics, depth summaries, and optional variant call outputs using tools such as samtools and BCFtools [@Danecek:2021]. Minimum depth thresholds are enforced to support high-confidence consensus generation. ACI, Kraken2, and IGV-Reports are also included to provide insight into amplicon metrics, contamination, and variant visualizations [@Young:ACI; @Wood:2019; @igv-reports:2024].
* Downstream Interpretation: Executes pathogen-specific analyses, including lineage and clade assignment with Pangolin and Nextclade [@Aksamentov:2021; @OToole:2021], wastewater lineage deconvolution with Freyja [@Karthikeyan:2022], and reference inspection with VADR [@Schaffer:2020].
* Phylogenetic Analysis (Optional): Supports multiple sequence alignment with MAFFT [@Katoh:2013] and phylogenetic tree construction with IQ-TREE [@Nguyen:2015] for comparative analyses and outbreak investigations. The pipeline utilizes snp-dists [@Seemann:snpdists] to generate SNP distance matrices, with heatcluster [@BeckstromSternberg:heatcluster] and phyTreeViz [@Moshi4:phytreeviz] integrated for the visual interpretation of clusters and phylogenetic relationships.

Each subworkflow is composed of modular Nextflow processes with containerized dependencies, allowing individual components to be enabled, disabled, or replaced through parameter configuration.

All results from individual Nextflow processes are summarized in a CSV file and MultiQC report [@Ewels:2016]. 

### Reproducibility and Portability

CECRET leverages containerization (Docker, Singularity [@Kurtzer:2017], or Apptainer) to ensure consistent execution across local systems, high-performance computing clusters, and cloud environments. All software dependencies are version-controlled within containers, minimizing variability due to system configuration.

The workflow includes predefined execution profiles (e.g., local, Docker, Slurm, test datasets) and parameter schema validation to enforce correct usage and improve user experience. Continuous integration testing is implemented to validate pipeline behavior across multiple configurations.

### Operation in Restricted and Offline Environments

Public health laboratories often operate in environments with limited or restricted internet access. To support these constraints, CECRET provides mechanisms to decouple runtime execution from external data dependencies:

* Users may supply predownloaded Nextclade datasets to avoid runtime downloads.
* Containerized tools (e.g., Pangolin, Kraken2, Freyja) include bundled or version-pinned resources to ensure consistent lineage assignment.
* Optional parameters allow disabling of database updates to maintain reproducibility in controlled environments.

These features enable reliable execution in secure or air-gapped systems while preserving analytical consistency.

The reliance on cloud-based assets and remote plugins in some modern workflows can introduce challenges in restricted-access or high-security laboratory environments. CECRET was designed with operational portability in mind, minimizing external dependencies and ensuring an execution model that is particularly suited for air-gapped or network-restricted public health workstations.

## Research Impact Statement

CECRET has been used as a production workflow for viral genomic surveillance in public health settings, supporting routine analysis of amplicon sequencing data during the SARS-CoV-2 pandemic and subsequent outbreaks.

At the Utah Public Health Laboratory (UPHL), the pipeline was deployed to generate consensus genomes and lineage assignments from Illumina-based amplicon sequencing, contributing to state-level surveillance efforts (BioProject PRJNA614995). The workflow has since been adapted to additional pathogens, including MPOX (PRJNA843095) and Measles Virus (PRJNA1293457), demonstrating its applicability to emerging public health threats.

The software is distributed as part of the StaPH-B toolkit [@Florek:staphb], facilitating adoption across public health laboratories and providing integration within a broader ecosystem of bioinformatics tools. The repository includes predefined test datasets, multiple execution profiles, and continuous integration workflows that support reproducibility and validation across environments.

CECRET is actively maintained, with regular updates to ensure compatibility with evolving lineage classification tools and to incorporate updated primer schemes and pathogen-specific configurations. The pipeline’s modular design has enabled rapid adaptation to new surveillance contexts, including wastewater-based epidemiology through integration with lineage deconvolution tools.

Together, these features support the use of CECRET as a stable, reproducible, and operationally focused workflow for high-throughput viral genomic surveillance.

## Validation and Benchmarking

Consensus sequences for 21 amplicon sequencing samples (15 SARS-CoV-2, 6 Measles) generated by CECRET (version 3.72.26090) and nf-core/viralrecon (version 2.6.0) were aligned to a common reference genome to ensure positional consistency. Pairwise SNP differences were calculated by counting positions at which two sequences differed among unambiguous nucleotides (A, C, G, T), excluding positions containing ambiguous bases (N) or gaps in either sequence. The proportion of ambiguous bases (% Ns) was calculated as the fraction of N characters relative to the full genome length. Lineage concordance was assessed using Pangolin and Nextclade, with agreement defined as identical lineage assignments between workflows.

### Table 1. Selected results for Samples Run with CECRET and nf-core/viralrecon
| Sample      | Category      |   Reads |   SNPs |   CECRET_N |   VR_N | Pango_C    | Pango_V    | Next_C      | Next_V      |
|:------------|:--------------|--------:|-------:|-----------:|-------:|:-----------|:-----------|:------------|:------------|
| SRR38049469 | Measles       |  621680 |      1 |       1.74 |   1.35 | -          | -          | D8          | D8          |
| SRR38049468 | Measles       |  545404 |      1 |      11.66 |   6.07 | -          | -          | D8          | D8          |
| SRR38049458 | Measles       |  765014 |      2 |       1.2  |   0.39 | -          | -          | D8          | D8          |
| SRR38049457 | Measles       |  692631 |      1 |       1.11 |   0.93 | -          | -          | D8          | D8          |
| SRR38049456 | Measles       |  547665 |      2 |       0.61 |   0.39 | -          | -          | D8          | D8          |
| SRR37335748 | Measles       |  424030 |      1 |       9.71 |   4.64 | -          | -          | B3          | B3          |
| SRR38006100 | High coverage |  299159 |      7 |      22.64 |   4.22 | XFG.2.5    | XFG.2.5.1  | 25C         | 25C         |
| SRR38006099 | High coverage |  537532 |      0 |       1.34 |   0.56 | XFG.23.1.4 | XFG.23.1.4 | 25C         | 25C         |
| SRR38006098 | High coverage | 1659941 |      1 |       1.68 |   0.35 | XFG.14.1   | XFG.14.1   | 25C         | 25C         |
| SRR38006097 | High coverage | 1420309 |      0 |       1.54 |   0.33 | XFG.2.5.1  | XFG.2.5.1  | 25C         | 25C         |
| SRR38006096 | High coverage |  873801 |      0 |       7.73 |   1.62 | XFG.6.2.2  | XFG.6.2.2  | 25C         | 25C         |
| SRR38006095 | High coverage |  233207 |      0 |      31.9  |   9.47 | Unassigned | XFG.3.6    | 25C         | 25C         |
| SRR38006093 | High coverage |  174702 |      0 |      30.68 |  10.47 | Unassigned | XFZ.2      | recombinant | recombinant |
| SRR38006092 | High coverage |  318099 |      0 |       5.49 |   1.43 | XFG.3.6    | XFG.3.6    | 25C         | 25C         |
| SRR38006091 | High coverage |  695713 |      0 |      10.03 |   3.15 | XFG.1.1    | XFG.1.1    | 25C         | 25C         |
| SRR38006090 | High coverage |  817668 |      0 |       0.73 |   0.6  | PQ.2.1.6   | PQ.2.1.6   | 25B         | 25B         |
| SRR36104323 | Low coverage  |   12482 |    nan |      99.5  | 100    | Unassigned | N/A        | N/A         | N/A         |
| SRR36104328 | Low coverage  |  102505 |    nan |     100    | 100    | Unassigned | N/A        | N/A         | N/A         |
| SRR36104302 | Low coverage  |  461127 |      0 |      94.68 |  88.59 | Unassigned | Unassigned | N/A         | 25C         |
| SRR36020460 | Low coverage  |  671807 |      0 |      96.07 |  93.08 | Unassigned | Unassigned | N/A         | N/A         |
| SRR35862479 | Low coverage  |  750481 |      0 |      97.83 |  97.83 | Unassigned | Unassigned | N/A         | N/A         |

### Table 2. Performance by Sample Type
| Category      |   N |   Med Reads |   Med SNP Diff |   Med % N (CECRET) |   Med % N (VR) | Concordance   |
|:--------------|----:|------------:|---------------:|-------------------:|---------------:|:--------------|
| High coverage |  10 |      616622 |              0 |               6.61 |           1.52 | 70.0%         |
| Low coverage  |   5 |      461127 |              0 |              97.83 |          97.83 | 100.0%        |
| Measles       |   6 |      584672 |              1 |               1.47 |           1.14 | 100.0%        |

CECRET was evaluated on 21 amplicon sequencing samples (15 SARS-CoV-2, 6 Measles) spanning a range of coverage profiles. Consensus sequences were compared against those generated using nf-core/viralrecon (Table 1). Across all 19 samples that successfully completed both pipelines (Table 2), consensus sequences differed by a median of 0 SNPs (range: 0–7). CECRET successfully generated consensus for all 21 samples, while nf-core/viralrecon failed on two low-coverage samples (<100k reads). Pangolin lineage concordance was 76.9%, with mismatches primarily due to "Unassigned" calls in high-N samples, while Measles clade concordance was 100%.

## AI Usage Disclosure

Generative AI tools, specifically ChatGPT and Google Gemini, were utilized during various stages of the development, maintenance, and documentation of the CECRET workflow.
The specific applications of AI included:
* Troubleshooting and Debugging: AI models were used to identify and resolve complex runtime errors and issues encountered during the pipeline's iterative development.
* Architectural Refinement: AI assisted in refining the workflow architecture to maintain compatibility with evolving Nextflow standards. This included major refactoring such as the migration across different Nextflow versions and the removal of deprecated environment output channels to improve workflow efficiency and structure.
* nf-core Compatibility: AI assisted in ensuring the workflow adhered to nf-core standards, specifically helping to resolve issues related to schema validation and the integration of linting tools into the continuous integration (CI) suite.
* Documentation and User Interface: The help text provided within the workflow’s command-line interface was generated and refined using AI to improve clarity for the end user.
* Manuscript Preparation: AI was used to generate the initial rough outline of this paper to ensure all necessary scholarly components were addressed.

Verification and Quality Control: All code adjustments, architectural changes, and documentation text suggested by AI were manually reviewed, tested, and verified for technical accuracy and biological relevance by the maintainer. The final consensus generation logic and tool integration remains grounded in established bioinformatics best practices.
