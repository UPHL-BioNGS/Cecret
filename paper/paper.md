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
 - name: Huntsman Cancer Institute, Department of Oncological Sciences or Medicine?, University of Utah, Salt Lake City, Utah, USA
   index: 2
date: 24 April 2026
bibliography: paper.bib
---

## Summary

CECRET is an open-source bioinformatics workflow for reference-based consensus genome generation from amplicon sequencing data. The pipeline was developed at the Utah Public Health Laboratory (UPHL) to address a specific operational need: producing high-quality consensus sequences from Illumina data generated using ARTIC-style amplicon protocols.

While widely adopted bioinformatics workflows exist, many were initially optimized for Oxford Nanopore Technologies (ONT) data or designed as general-purpose viral analysis pipelines. In contrast, CECRET is purpose-built for amplicon-based Illumina sequencing, where accurate primer trimming, depth-aware variant calling, and standardized downstream lineage assignment are critical for reliable consensus generation.

The workflow automates the transformation of raw sequencing reads into analysis-ready outputs through a modular pipeline that includes read preprocessing, reference alignment, primer trimming, consensus generation, and downstream pathogen typing. It integrates established tools such as BWA or minimap2 for alignment, iVar for primer trimming and consensus generation, and lineage/clade assignment tools including Nextclade and Pangolin. Optional modules support multiple sequence alignment and phylogenetic inference for outbreak investigation.

CECRET is implemented in Nextflow and distributed with containerized dependencies (Docker, Singularity, or Apptainer), enabling reproducible execution across local, high-performance computing, and cloud environments. Although initially developed for SARS-CoV-2, the workflow has been extended to support additional viral pathogens such as mpox and measles, and can be adapted to other small viral genomes given appropriate reference and primer scheme inputs.


## Statement of Need

Rapid and reproducible genomic surveillance is essential for monitoring viral evolution and informing public health responses during outbreaks [source]. Amplicon-based sequencing approaches, particularly those derived from ARTIC protocols, have been widely adopted due to their scalability and cost efficiency. However, the bioinformatic processing of these data presents distinct challenges that are not consistently addressed by general-purpose workflows.

At the onset of the SARS-CoV-2 pandemic, existing bioinformatics pipelines provided by the ARTIC Network were primarily optimized for Oxford Nanopore Technologies (ONT) data [source]. Public health laboratories frequently implemented hybrid workflows combining ARTIC primer schemes with Illumina sequencing platforms [source], creating a gap in available tooling. In these settings, accurate consensus generation requires post-alignment primer trimming and depth-aware variant filtering to avoid technical artifacts introduced during amplification.

CECRET was developed to address this gap by providing a standardized, reference-based workflow tailored specifically to amplicon-based Illumina data. The pipeline encodes best practices for:

* post-alignment primer trimming to remove amplification artifacts,
* minimum depth thresholds for high-confidence consensus generation,
* and automated integration of lineage and clade assignment tools.

Unlike broader viral analysis workflows (e.g., nf-core/viralrecon), which support diverse sequencing strategies including shotgun metagenomics and *de novo* assembly, CECRET prioritizes a constrained, reference-guided approach optimized for high-throughput public health surveillance. This design enables consistent, reproducible outputs with minimal parameter tuning, making it suitable for routine operational use.

The primary users of CECRET are public health laboratories, clinical genomics groups, and researchers conducting targeted viral surveillance using amplicon sequencing. The workflow is particularly suited to scenarios requiring rapid turnaround and standardized outputs, such as outbreak response and longitudinal monitoring programs.


## State of the Field

Early in the SARS-CoV-2 pandemic, viral genomic surveillance workflows were often assembled from heterogeneous tools with limited standardization. While the ARTIC Network provided widely adopted laboratory protocols and associated bioinformatics guidance, these pipelines were primarily designed for ONT data and required adaptation for other sequencing platforms [source]. Public health laboratories using Illumina-based amplicon sequencing frequently relied on custom scripts and locally maintained workflows, resulting in variability in analytical outcomes and reduced reproducibility across institutions.

Since 2020, the field has progressed toward standardized, containerized workflows that emphasize reproducibility and portability. Frameworks such as Nextflow and nf-core have enabled the development of community-maintained pipelines that support a wide range of sequencing modalities and analytical goals [source]. These include comprehensive workflows capable of handling metagenomic data, reference-based analyses, and *de novo* assembly within a single framework.

Despite these advances, there remains a need for workflows that prioritize operational simplicity and encode domain-specific best practices for targeted surveillance applications. In particular, amplicon-based Illumina sequencing introduces methodological constraints—such as the requirement for post-alignment primer trimming—that are not consistently enforced in general-purpose pipelines.

CECRET addresses this niche by providing a focused, reference-based workflow optimized for small viral genomes and high-throughput surveillance contexts. Rather than maximizing flexibility across sequencing strategies, the pipeline emphasizes robustness, reproducibility, and minimal configuration for a well-defined class of analyses commonly performed in public health laboratories.


## Build vs. Contribute Justification

CECRET was developed as a standalone workflow rather than as a contribution to an existing pipeline due to the specific operational constraints of public health genomic surveillance during the early stages of the SARS-CoV-2 pandemic.

At that time, available workflows were either narrowly tailored to ONT data (e.g., ARTIC pipelines) or designed as flexible research frameworks supporting multiple sequencing strategies. These approaches did not provide a stable, Illumina-focused solution with enforced best practices for amplicon data processing. In particular, critical steps such as post-alignment primer trimming and depth-based consensus filtering were not consistently implemented as default behaviors.

Developing a dedicated workflow enabled the encoding of these domain-specific requirements as fixed or strongly guided defaults, reducing the need for local customization and minimizing variability between laboratories. This design approach prioritizes reproducibility and consistency over configurability, which is a key requirement in public health settings where standardized outputs are necessary for downstream reporting and data sharing.

In addition, CECRET integrates pathogen-specific downstream tools (e.g., lineage and clade assignment, wastewater deconvolution) into a single pipeline, reducing the need for manual chaining of independent tools. This level of integration would have required substantial restructuring of existing general-purpose workflows.

As the field has matured, comprehensive pipelines such as nf-core/viralrecon have expanded in scope and capability. However, CECRET continues to serve a complementary role by providing a constrained, production-oriented workflow optimized for high-throughput amplicon-based surveillance using Illumina data.


## Software Design

CECRET is implemented using the Nextflow workflow manager (DSL2), enabling modular composition, reproducibility, and portability across diverse computational environments. The pipeline is organized into discrete subworkflows and processes that reflect the logical stages of amplicon-based consensus generation.

### Architecture and Modularity

The workflow is structured into several primary subworkflows:

* **Initialization**: Standardizes input formats (paired-end Illumina reads, single-end reads, nanopore reads, and FASTA inputs), validates parameters, and prepares reference genome resources.
* **Consensus Generation**: Performs read preprocessing, alignment, and primer trimming. Illumina reads are aligned using BWA by default (with optional minimap2 support), followed by primer trimming and consensus generation using iVar. For nanopore data, ARTIC-based tools are used.
* **Quality Control and Variant Calling**: Generates coverage metrics, depth summaries, and optional variant call outputs using tools such as samtools and BCFtools. Minimum depth thresholds are enforced to support high-confidence consensus generation.
* **Downstream Interpretation**: Executes pathogen-specific analyses, including lineage and clade assignment (e.g., Pangolin, Nextclade) and, where applicable, wastewater lineage deconvolution (Freyja).
* **Phylogenetic Analysis (Optional)**: Supports multiple sequence alignment (MAFFT) and phylogenetic tree construction (IQ-TREE) for comparative analyses and outbreak investigations.

Each subworkflow is composed of modular Nextflow processes with containerized dependencies, allowing individual components to be enabled, disabled, or replaced through parameter configuration.

### Reproducibility and Portability

CECRET leverages containerization (Docker, Singularity, or Apptainer) to ensure consistent execution across local systems, high-performance computing clusters, and cloud environments. All software dependencies are version-controlled within containers, minimizing variability due to system configuration.

The workflow includes predefined execution profiles (e.g., local, Docker, Slurm, test datasets) and parameter schema validation to enforce correct usage and improve user experience. Continuous integration testing is implemented to validate pipeline behavior across multiple configurations.

### Operation in Restricted and Offline Environments

Public health laboratories often operate in environments with limited or restricted internet access. To support these constraints, CECRET provides mechanisms to decouple runtime execution from external data dependencies:

* Users may supply predownloaded Nextclade datasets to avoid runtime downloads.
* Containerized tools (e.g., Pangolin, Freyja) include bundled or version-pinned resources to ensure consistent lineage assignment.
* Optional parameters allow disabling of database updates to maintain reproducibility in controlled environments.

These features enable reliable execution in secure or air-gapped systems while preserving analytical consistency.

The reliance on cloud-based assets and remote plugins in some modern workflows can introduce challenges in restricted-access or high-security laboratory environments. CECRET was designed with operational portability in mind, minimizing external dependencies and ensuring a 'self-contained' execution model that is particularly suited for air-gapped or network-restricted public health workstations.

## Research Impact Statement

CECRET has been used as a production workflow for viral genomic surveillance in public health settings, supporting routine analysis of amplicon sequencing data during the SARS-CoV-2 pandemic and subsequent outbreaks.

At the Utah Public Health Laboratory (UPHL), the pipeline was deployed to generate consensus genomes and lineage assignments from Illumina-based amplicon sequencing, contributing to state-level surveillance efforts. The workflow has since been adapted to additional pathogens, including mpox, demonstrating its applicability to emerging public health threats.

The software is distributed as part of the StaPH-B toolkit [source], facilitating adoption across public health laboratories and providing integration within a broader ecosystem of bioinformatics tools. The repository includes predefined test datasets, multiple execution profiles, and continuous integration workflows that support reproducibility and validation across environments.

CECRET is actively maintained, with regular updates to ensure compatibility with evolving lineage classification tools (e.g., Pangolin, Nextclade) and to incorporate updated primer schemes and pathogen-specific configurations. The pipeline’s modular design has enabled rapid adaptation to new surveillance contexts, including wastewater-based epidemiology through integration with lineage deconvolution tools.

Together, these features support the use of CECRET as a stable, reproducible, and operationally focused workflow for high-throughput viral genomic surveillance.

## Validation and Benchmarking

TO BE ADDED!!!
CECRET was evaluated on 24 SARS-CoV-2 amplicon sequencing samples spanning a range of coverage profiles. Consensus sequences were compared against those generated using nf-core/viralrecon and the artic Illumina workflow. Across all samples, lineage assignments were identical, and consensus sequences differed by a median of 1 SNP (range: 0–3). The proportion of ambiguous bases was comparable between workflows, with a median of 2.1% for CECRET and 2.3% for viralrecon. These results demonstrate that CECRET produces concordant outputs while enforcing standardized processing steps for amplicon-based Illumina data.

Consensus sequences generated by CECRET, nf-core/viralrecon, and the ARTIC workflow were aligned to a common reference genome to ensure positional consistency. Pairwise SNP differences were calculated by counting positions at which two sequences differed among unambiguous nucleotides (A, C, G, T), excluding positions containing ambiguous bases (N) or gaps in either sequence. The proportion of ambiguous bases (% Ns) was calculated as the fraction of N characters relative to the full genome length. Lineage concordance was assessed using Pangolin, with agreement defined as identical lineage assignments between workflows.

Operational benchmarking highlighted differences in input schema resilience. We observed that strict validation layers in specialized pipelines—while beneficial for ensuring data quality in standardized runs—can create execution barriers when input formats deviate from expected templates (e.g., auxiliary columns in primer BED files). CECRET prioritizes interoperability by utilizing industry-standard tool behaviors that maintain performance across a wider variety of input formatting variations.

I think I just need ~25 samples (20 SARS-CoV-2 where 10 are good, 5 are mediocre, 5 are bad, and 5 are measles). Metrics are % identity to each other, SNP differences, mean depth or % coverage and pangolin and nextclade assignments.

To evaluate performance, CECRET was compared against an established reference-based workflow using a set of amplicon sequencing samples spanning a range of coverage profiles. Summary metrics are shown in Table 1, with performance stratified by sample type in Table 2.

### Table 1.
| Metric                     | CECRET | nf-core/viralrecon | artic |
| -------------------------- | ------ | - | - |
| Number of samples          |        |  |  |
| Median SNP difference      | —      |   |  |
| SNP difference (range)     | —      |  |  |
| Lineage concordance (%)    |        |   |  |
| Median % Ns                |        |  |  |
| % Ns (range)               |        |  |  |
| Median genome coverage (%) |        |  |  |

### Table 2.
| Sample Type      | N | Median SNP Diff (vs viralrecon) | Median SNP Diff (vs ARTIC) | Median % Ns (CECRET) | Median % Ns (viralrecon) | Median % Ns (ARTIC) | Lineage Concordance (vs viralrecon) | Lineage Concordance (vs ARTIC) |
| ---------------- | - | ------------------------------- | -------------------------- | -------------------- | ------------------------ | ------------------- | ----------------------------------- | ------------------------------ |
| High coverage    |   |                                 |                            |                      |                          |                     |                                     |                                |
| Low coverage     |   |                                 |                            |                      |                          |                     |                                     |                                |
| Amplicon dropout |   |                                 |                            |                      |                          |                     |                                     |                                |



## AI Usage Disclosure

Generative AI tools, specifically ChatGPT and Google Gemini, were utilized during various stages of the development, maintenance, and documentation of the CECRET workflow.
The specific applications of AI included:
* Troubleshooting and Debugging: AI models were used to identify and resolve complex runtime errors and issues encountered during the pipeline's iterative development.
* Architectural Refinement: AI assisted in adjusting the workflow's architecture to maintain compatibility with evolving Nextflow standards. This included major refactoring such as the migration across different Nextflow versions and the removal of deprecated env output channels to improve workflow efficiency and structure.
* nf-core Compatibility: AI assisted in ensuring the workflow adhered to nf-core standards, specifically helping to resolve issues related to schema validation and the integration of linting tools into the continuous integration (CI) suite.
* Documentation and User Interface: The help text provided within the workflow’s command-line interface was generated and refined using AI to improve clarity for the end user.
* Manuscript Preparation: AI was used to generate the initial rough outline of this paper to ensure all necessary scholarly components were addressed.

Verification and Quality Control: All code adjustments, architectural changes, and documentation text suggested by AI were manually reviewed, tested, and verified for technical accuracy and biological relevance by the maintainer. The final consensus generation logic and tool integration remains grounded in established bioinformatics best practices.
