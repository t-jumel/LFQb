# LFQ benchmark
Benchmarking label-free quantification (LFQ) in bottom-up proteomics by DIA-LC-MS and DIA-NN.


This script analyzes DIA-NN precursor and protein group matrices of LFQ benchmarks resulting
in plots and summary statistics, which provide a meaningful comparison of multiple result sets 
(e.g. different DIA-NN settings) or troubleshooting utlity.
   
    
      


This script analyzes LFQ benchmark datasets and serves as alternative to scripts and packages
from the following publications:

- Kuharev, Jörg, et al. "In‐depth evaluation of software tools for data‐independent acquisition based label‐free quantification." Proteomics 15.18 (2015): 3140-3151.

- Navarro, Pedro, et al. "A multicenter study benchmarks software tools for label-free proteome quantification." Nature biotechnology 34.11 (2016): 1130-1136.



A nice resource for benchmark raw files from various instrument types can be found here:
- https://www.ebi.ac.uk/pride/archive/projects/PXD028735

# Quick Start Guide
## Samples

<img src="readme_figures/01.png" alt="sample mixtures" width="250"/>

- Sample mixtures tpyically derived from commercial digests.
- Aim at ca. 3 replicates per condition.
- This script should handle different expected fold-changes and even 2-species mixtures directly, the layout (yeast upregulated from B to A etc.) is hard-coded, deviating might result in some summary stats to result in nonsense values.
- Raw data processing with DIA-NN against SwissProt fasta is recommended and only the respective matrices can be directly analysed.


## Analysis Preparation

### R packages
Install packages, CRAN automatically, Bioconductor below.
```
if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(
  tidyr,
  matrixStats,
  mratios,
  magrittr,
  statmod,
  scrime,
  moments,
  reshape2,
  grid,
  gridExtra,
  ggplot2,
  cowplot,
  scales)
```
```
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("Biobase")
```
### Input Requirements
Prepare a folder containing exactly one precursor and one protein group matrix in .tsv format
with "pr_matrix" and "pg_matrix" in the respective filename. Other filetypes are ignored. 
(Basically the results of one DIA-NN search, but additional "first-pass" files need to be moved away)

<img src="readme_figures/02.png" alt="sample mixtures" width="300"/>

### Variables and Filter settings
In the variables section, set appropriate parameters. The most important ones can be seen here,
such as the folder location, column indices of the LFQ values, and filter settings.
The example below corresponds to:
- precursor and protein groups will be counted as "IDs" if a value is reported in at least 2 of 4 replicates in both condition A and B
- precursor and protein groups will be counted as "Quantified" if they are identified and have a CV less than 20% in both condition A and B
- quantified protein groups are subjected to differential expression analysis with limma with a statistical cutoff of 1% (0.01),
 however, protein groups are only recognized as up- or downregulated if the log2 fold-change exceeds +-0.5

```
folder_input <- "C:/Users/Tobias/Desktop/Test_Input"

cond_ctr <- "LFQ_B"
cond_exp <- "LFQ_A"

# # Pretyped 4 replicates per sample.
col_exp_Prot <- 6:9
col_ctr_Prot <- 10:13
col_exp_Prec <- 11:14
col_ctr_Prec <- 15:18

# Filter variables, listed in output folder name.
limit_MV <- (2 / 4)
limit_CV <- 20
limit_FC <- 0.5

# p_adj cut-off for diff. expr. analysis by limma.
alpha_limma <- 0.01
```

