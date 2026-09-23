# Patient-Specific Molecular Network Entropy Defines Two Prognostic Ecosystem States in Lung Adenocarcinoma

This repository contains all scripts used to generate the analyses, figures, tables, and results reported in the manuscript:

**"Patient-Specific Molecular Network Entropy Defines Two Prognostic Ecosystem States in Lung Adenocarcinoma."**

The workflow reconstructs patient-specific gene co-expression networks from TCGA-LUAD using the LIONESS framework, quantifies multiscale network entropy across curated biological modules, identifies ecosystem states through unsupervised clustering, and evaluates their molecular, immune, network, and prognostic characteristics. A frozen classifier derived from TCGA was subsequently applied to the independent GSE72094 cohort for external validation without reclustering or model retraining.

## Data Sources

- TCGA-LUAD (discovery cohort; n = 517)
- GSE72094 (external validation cohort)

## Main Analyses

- Data preprocessing and quality control
- Patient-specific network construction using LIONESS
- Module-level network entropy calculation
- Ecosystem-state discovery and stability assessment
- Pathway and immune microenvironment characterization
- Mutation and network topology analyses
- Survival modeling and external validation

## Author

**Bikram Sahoo, Ph.D.**  
Email: biks.kestro@gmail.com
