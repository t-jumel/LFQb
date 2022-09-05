# Overview ------------------------------------------------------------

### Literature:
# Demichev, Vadim, et al. "DIA-NN: neural networks and interference correction enable deep proteome coverage in high throughput." Nature methods 17.1 (2020): 41-44.
# Kuharev, Jörg, et al. "In-depth evaluation of software tools for data-independent acquisition based label-free quantification." Proteomics 15.18 (2015): 3140-3151.
# Navarro, Pedro, et al. "A multicenter study benchmarks software tools for label-free proteome quantification." Nature biotechnology 34.11 (2016): 1130-1136.
# Ritchie, Matthew E., et al. "limma powers differential expression analyses for RNA-sequencing and microarray studies." Nucleic acids research 43.7 (2015): e47-e47.



### Disclaimer:
# Hi, my name is tobias and I hope you can make use
# of this script. Please excuse if the code is not elegant
# as I am a beginner with R. I am thankful for suggestions.
# This script was inspired by the "classical LFQbench papers"
# (Kuharev et al. and Navarro et al.), 
# and serves as an alternative for the related packages and scripts.



### Contact:
# Tobias Jumel
# jumel@mpi-cbg.de
# Shevchenko Lab, MPI-CBG
# 05.09.2022
# Limma contributed by Andre Gohr, MPI-CBG scientific computing.



### License: (GPL-3)
# <LFQbenchmark script to process DIA-NN output for label-free
# quantification benchmarks in bottom-up proteomics by mass spectrometry >
# Copyright (C) <2022> <Tobias Jumel>
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public
# License along with this program. If not, see
# <http://www.gnu.org/licenses/>





# To run this Script ------------------------------------------------------

# 1) Check package installations, especially bioconductor ones.


# 2) Select a folder location (as variable "folder_input").
# I can recommend "path copy copy".
# This folder must contain exactly one protein group
# and precursor .tsv matrix with DIA-NN style column names
# with "pg_matrix" and "pr_matrix" in the respective filenames.
# Usage of +MBR requires the undesired first-pass or 
# second pass matrices to be moved away , e.g. into a subfolder.


# 3) Adjust all parameters under "variables" section, in most cases
# this applies just to folder_input, filter settings, and column indices.
# Variables of species not contained in the data are ignored.

# This script assumes human to be stable, yeast up, and both E.coli
# and C.elegans to be ownregulated. Any deviation requires minor adaption of
# the function f.LFQbenchmark.limma.interpretation according 
# to the examples, otherwise the confusion matrix summary stats become nonsense.



# 4) For paired stat. test in f.limma.p_adj just enable "reps" and "limma.rep" (3x enable)
# By default neither are used for an unpaired comparison of conditions.


# Then execute and find the output folder (csv files, plots)
# within the specified input folder.



# Introduction -------------------------------------------------------------

### General
# This script analyses benchmarks of peptide and protein group
# label-free quantification by bottom-up proteomics with DIA LC-MS and DIA-NN.



### Samples
# The benchmarks use 2 samples, each consisting of multiple species
# mixed together in different percentages.(see Literature)
# classical Sample A: 65% Human, 30% Yeast, 5% E.coli
# classical Sample B: 65% Human, 15% Yeast, 20% E.coli
# Often measured with ca. 3 technical replicates.

# Note: This script works directly on 2-species (Human,Yeast) and 
# 4-species benchmarks (Human, Yeast, E.coli, C.elegans).
# Please note instructions in the "To run this script" section



### Analysis
# When comparing peptide and protein group fold-changes from sample B to A,
# one expects Human entries to be stable, Yeast 2x up, and E.coli
# 4x down-regulated (non-log scale, based on sample mixtures). The benchmarks reveal to what extent 
# measured log2 fold-changes deviate from expected log2 fold-changes,
# including e.g. global over-or underestimation and identification errors
# that often remain invisible by not having expectation values available.
# This script has a strong focus on visualizing and measuring the degree
# of log2 fold-change over-or underestimation, (also called ratio compression 
# or ratio expansion), mainly with the Asymmetry_Factor.


# Parts of this script are hard-coded matching the given example.
# (function f.LFQbenchmark.limma.interpretation)
# While expected fold-changes can be altered under variables to match different samples,
# other parts of the script are hard-coded expecting human protein groups
# to be stable, yeast ones to be up-regulated and E.coli ones to be down-regulated
# from Sample B to Sample A. Diverting from this orientation might result
# in some summary statistics to become nonsense values.



### Limma and LFQ benchmarks

# Apart from just matching measured against expected log2 fold-changes
# the protein groups are subjected to differential expression analysis
# using limma. It is measured which protein groups are up-regulated,
# down-regulated, or stable. These results are matched against the 
# expectations (all Human ones are stable, all Yeast ones are up-regulated,
# all E.coli ones are down-regulated).

# Examples: 
# A Yeast protein group measured as up-regulated is a true positive.
# A Human protein group measured as stable is a true-negative.
# A human protein group measured as up-regulated is a false positive.
# A Yeast protein group measured as stable is a false negative.

# The count of true positives etc. equals a "confusion matrix" layout.
# The true positive count and the false positive rate (deFDR)
# are the most important summary statistics to describe the performance
# and reliability of a workflow. They allow for better optimization of
# e.g. DIA-NN settings rather than just maximizing the raw proteome coverage.



### Interpretation

# 1) Check all plots for anything unusual incl Precursor-level. Expect the unexpected.


# 2) Check for appropriate normalisation.
# Potential erroneous offsets (medians or scatter plots) could invalidate the
# summary stats generated. While ratio compression or expansion can move 
# Yeast and E.coli medians, the critical normalisation related errors typically shift 
# Human and Yeast medians or data points by the same distance and direction 
# on the log2 fold-change axis.


# 3) Check reliability to prevent excess false positives in differential expression analysis.
# Reliability is related to the deFDR value and the presence of aberrant quantifications
# in facet plots. Artificially increased protein group FDR leads to clear clusters
# of human protein groups around the expected Yeast log2 fold-change, and vice versa.
# A deFDR below 1 %, the absence of such clusters and no strong ratio expansion
# speak for excellent reliability. A workflow optimized to minimize these error sources
# might be less likely to result in false positives during differential expression applications.
# Depending on the scope and replicate numbers, deFDR values higher than 1 %
# might be tolerated.


# 4) Check fold-change over-or underestimation.
# One of the most overlooked issues in quantitative proteomics.
# Orbitrap data might result in ratio expansion while ToF data might result in the 
# opposing ratio compression, both cases depend on the analysis settings.
# A separate data frame and output file lists multiple stats related to these.
# The most useful one seems to be the Asymmetry_Factor. 

# A Factor of 1 means perfect symmetry and no bias. 
# 0.5 indicates a strong and undesirable underestimation, leading to slightly differentially
# abundant proteins not recognized correctly as such.
# The increased number of false negatives means
# differential expression results can be reliable but less sensitive.
# An Asymmetry_Factor of 2 indicates a strong and undesired fold-change overestimation.
# Might lead to a slightly increased deFDR and false positive rate 
# in differential expression applications.


# 5) Check Precision.
# Good Precision is a requirement for overall accuracy.
# When working with ca. 3 replicates and a CV filter of 20%,
# it seems well sufficient if a result set ends up with average and median
# CV at or below 5 %. The value aimed at might vary based on replicate number,
# proteome coverage, and intended application outside of benchmarks.
# While workflows with higher values might be reliable, a loss of 
# sensitivity in differential expression applications is expected at some point.


# 6) When optimizing, DIA-NN or instrument settings, etc.,
# aim to stay within a specified deFDR and ratio compression / extension
# while maximizing the true positive count. A TP count exceeding 2500 is quite good, 
# 3000 seems hard to reach and might be close to the upper limit. The scaling is non-linear, 
# increasing from 2.5k to 3k is more difficult than from 2k to 2.5k.





# Environment preparation -------------------------------------------------------

# Clear environment, RStudio plots, console, and packages.
# Will not clear package namespaces imported by other packages.
rm(list = ls())
dev.off(dev.list()["RStudioGD"])
invisible(lapply(
  paste0('package:', names(sessionInfo()$otherPkgs)),
  detach,
  character.only = TRUE,
  unload = TRUE
))
cat("\014")



# Install and load required base packages automatically.
# Bioconductor ones are separately below.
if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(
  tidyr,
  dplyr,
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



# Have updated bioconductor packages.
# Update or install often prohibited by RStudio not allowed to write in folders.
# Try open RStudio "run as administrator" and use .libPaths() to find and manually
# delete folders via windows file explorer to allow fresh re-install. 


### Pre-typed to install bioconductor packages.
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("limma")
# BiocManager::install("Biobase")

library(limma)
library(Biobase)




# Variables ---------------------------------------------------------------

# Explanation below at section end for better spacing and overview.


# Set working directory (just pre-typed, might not be used).
# setwd()

# Script version and name used in output filenames.
script_name <- "LFQb"
script_version <- "v3.2"


# Recommend use of "path copy copy".
# See "To run this Script" section for input requirements.
folder_input <- "//fileserver/project_jumel/All_LFQ_benchmarks/Reanalysis_external/Reanalysis_LFQbench_PXD028735/diaPASEF/diaPASEF06"



# Naming the 2 different benchmark sample mixtures.
# This is just for record-keeping and preventing manual errors 
# when assigning columns to input categories below.
cond_ctr <- "LFQ_B"
cond_exp <- "LFQ_A"


# From input tables, provide indices of columns
# with quantitative values for the 2 conditions named above.

# # Pretyped 3 replicates per sample type.
col_exp_Prot <- 6:8
col_ctr_Prot <- 9:11
col_exp_Prec <- 11:13
col_ctr_Prec <- 14:16

# # Pretyped 3 replicates from table of 4 replicates.
 # col_exp_Prot <- 6:8
 # col_ctr_Prot <- 10:12
 # col_exp_Prec <- 11:13
 # col_ctr_Prec <- 15:17

# # Pretyped 4 replicates per sample.
# col_exp_Prot <- 6:9
# col_ctr_Prot <- 10:13
# col_exp_Prec <- 11:14
# col_ctr_Prec <- 15:18


# Filter variables, listed in output folder name.
limit_MV <- (2 / 3)
limit_CV <- 20
limit_FC <- 0.5

rep_max <- (length(col_ctr_Prot) + length(col_exp_Prot)) / 2
rep_min <- rep_max * limit_MV


# p_adj cut-off for differential expression analysis by limma.
alpha_limma <- 0.01


# Expected fold-changes by sample mixtures.
expFC_human <- 0
expFC_yeast <- +1
expFC_ecoli <- -2
expFC_celegans <- -1

# Variables for plots.

# Colors
color_human <- "#199d76"
color_yeast <- "#d85f02"
color_ecoli <- "#7570b2"
color_celegans <-  "darkred"

# These plot limits are set a bit higher (0.2) than 
# needed to prevent value clipping with density plots.
# (Intentionally to prevent bug causing color loss)
FC_max <- +3.2
FC_min <- -3.2

# Export plot dimensions.
plot_height <- 10
plot_width <- 10
plot_res <- 900




### Explanation variables:
# script_name and script_version are used in filenames 
# of exported folders,tables and plots.

# folder_input is a folder containing one DIA-NN 
# protein group and precursor matrix to analyze.
# See "To run this script" section.

# cond_ctr and cond_exp contain names for the 2 conditions 
# (benchmark sample mixtures).

# col_exp_Prot etc. are the column indices for the 
# quantitative values in protein group and precursor matrices,
# in both cases for the experimental and control condition.

# limit_MV indicates the desired limit for missingness 
# of quantitative values. e.g. (2/4) or 0.5 means an entry 
# passes the respective filter step if at least
# 2 of 4 replicate measurements provide a non-0 and non-NA value
# in both conditions. Entries passing are called "IDs" or "identified".

# rep_max and rep_min indicate the range of value counts per condition.

# limit_CV indicates the limit of standard deviation in
# form of the coefficient of variation in %. 
# Entries are considered precise and pass the respective filter 
# if the CV is below e.g. 20 % in both conditions. 
# Entries passing are called "quantified".

# Measured log2 fold-changed between +- limit_FC are always 
# considered as "NOT sign. up- or down-regulated",
# no matter the limma result (p_adj).
# This value might be set to be between 0.5-1.0 .

# alpha_limma is the stat. significance cutoff 
# for the p_adj of limma. Always compare the resulting deFDR to this value,
# an excess false positive rate indicates potential errors in the data.

# expFC_ecoli and similar refer to the species-specific 
# expected log2 fold-changes.
# The values are based on mixing ratios of LFQbenchmark samples A and B.


# Variables for plots

# color_ecoli and similar are the species-specific colors for plots.

# FC_max and FC_min provide limits for plotting log2 fold-changes.

# Plot width and height in cm, resolution in dpi.





# Create output folder for plots and tables within input folder
# with filter variables in folder name. This allows running this script
# with various settings without overwriting prior results.
folder_output <-
  paste0(
    folder_input, "/",script_name,
    "_", script_version,
    "_", rep_min,"of", rep_max,
    "_CV", limit_CV,
    "_FC",limit_FC,
    "_a",alpha_limma)

dir.create(folder_output)



# Collect variables in data frame for export.
variables <- as.data.frame(
  c(
    script_name,
    script_version,
    folder_input,
    paste0("'", as.character(col_ctr_Prot), "'", collapse = ", "),
    paste0("'", as.character(col_exp_Prot), "'", collapse = ", "),
    paste0("'", as.character(col_ctr_Prec), "'", collapse = ", "),
    paste0("'", as.character(col_exp_Prec), "'", collapse = ", "),
    limit_MV,
    rep_max,
    rep_min,
    limit_CV,
    limit_FC,
    alpha_limma,
    cond_ctr,
    cond_exp,
    expFC_human,
    expFC_yeast,
    expFC_ecoli,
    expFC_celegans,
    color_human,
    color_yeast,
    color_ecoli,
    color_celegans,
    folder_output
  )
)

colnames(variables) <- ""
rownames(variables) <- c(
  "script_name",
  "script_version",
  "folder_input",
  "col_ctr_Prot",
  "col_exp_Prot",
  "col_ctr_Prec",
  "col_exp_Prec",
  "limit_MV",
  "rep_max",
  "rep_min",
  "limit_CV",
  "limit_FC",
  "alpha_limma",
  "cond_ctr",
  "cond_exp",
  "expFC_human",
  "expFC_yeast",
  "expFC_ecoli",
  "expFC_celegans",
  "color_human",
  "color_yeast",
  "color_ecoli",
  "color_celegans",
  "folder_output"
)



# Export variables as separate .csv file.
write.csv(variables,
          row.names = TRUE,
          paste0(folder_output, "/", "variables_",
                 basename(dirname(folder_output)),"_",
                 script_version,"_",
                 rep_min,"of", rep_max,
                 "_CV", limit_CV,
                 "_FC",limit_FC,
                 "_a",alpha_limma, ".csv"))



       
       


# Functions ---------------------------------------------------------------

# Recognize user defined functions (f.abc) and arguments (a.xyz)
# Plot functions are at the script end (f.p.scatter) together with
# the respective plot export function.



# Import DIA-NN matrix from specified folder into data frame.
# Uses partial filename match, e.g."pg_matrix", "pr_matrix".
# requires excess matrices to be moved away. (Applies when using MBR)
# e.g. Prot0 <- f.import.matrix(folder_input, "pg_matrix", "\t")
f.import.matrix <- function(a.folder,
                            a.matrix,
                            a.sep) {
  read.csv(sep = a.sep,
           check.names = FALSE,
           list.files(a.folder,
                      pattern = a.matrix,
                      full.names = TRUE))
}




# In DIA-NN matrices, rename columns containing quantitative values
# from filepath to short (stripped) replicate name.
# Remove filepath and suffix, remainder is used as replicate name.
# e.g. Prot1 <- f.strip.sample.column.name(Prot1, col_ctr_Prot)
f.strip.sample.column.name <- function(a.df, a.columns) {
  
  names(a.df)[a.columns] <- gsub(".raw$", "", names(a.df)[a.columns])
  names(a.df)[a.columns] <- gsub(".dia$", "", names(a.df)[a.columns])
  names(a.df)[a.columns] <- gsub(".mzML$", "", names(a.df)[a.columns])
  names(a.df)[a.columns] <- gsub(".wiff$", "", names(a.df)[a.columns])
  names(a.df)[a.columns] <- gsub(".wiff2$", "", names(a.df)[a.columns])
  names(a.df)[a.columns] <- gsub(".d$", "", names(a.df)[a.columns])
  
  names(a.df)[a.columns] <- gsub(".raw$", "", names(a.df)[a.columns])
  names(a.df)[a.columns] <- gsub(".dia$", "", names(a.df)[a.columns])
  names(a.df)[a.columns] <- gsub(".mzML$", "", names(a.df)[a.columns])
  names(a.df)[a.columns] <- gsub(".wiff$", "", names(a.df)[a.columns])
  names(a.df)[a.columns] <- gsub(".wiff2$", "", names(a.df)[a.columns])
  names(a.df)[a.columns] <- gsub(".d$", "", names(a.df)[a.columns])
  
  names(a.df)[a.columns] <- gsub(".raw$", "", names(a.df)[a.columns])
  names(a.df)[a.columns] <- gsub(".dia$", "", names(a.df)[a.columns])
  names(a.df)[a.columns] <- gsub(".mzML$", "", names(a.df)[a.columns])
  names(a.df)[a.columns] <- gsub(".wiff$", "", names(a.df)[a.columns])
  names(a.df)[a.columns] <- gsub(".wiff2$", "", names(a.df)[a.columns])
  names(a.df)[a.columns] <- gsub(".d$", "", names(a.df)[a.columns])
  
  names(a.df)[a.columns] <- basename(names(a.df)[a.columns])
  # names(a.df)[a.columns] <- sub(".[^.]+$", "", colnames(a.df)[a.columns])
  return(a.df)
}




# Remove entries of MaxQuant contaminant fasta from an e.g. DIA-NN pg_matrix.
# e.g. Prot1 <- f.remove.MQ.cont(Prot1)
f.remove.MQ.cont <- function(a.df) {
  a.df <- a.df[!grepl("TREMBL", a.df[, "Genes"]),]
  a.df <- a.df[!grepl("SWISS-PROT", a.df[, "Genes"]),]
  a.df <- a.df[!grepl("REFSEQ", a.df[, "Protein.Ids"]),]
  a.df <- a.df[!grepl("ENSEMBL", a.df[, "Protein.Ids"]),]
  a.df <- a.df[!grepl("H-INV:HIT", a.df[, "Protein.Ids"]),]
  a.df <- a.df[!grepl("Streptavidin", a.df[, "Protein.Ids"]),]
}




# For each entry and condition, calculate stats of the replicate measurements
# such as count of non-0 and non-NA values, mean, standard deviation, 
# coefficient of variation (CV - standard deviation in % of the mean)
# e.g. Prot1 <- f.row.stats(Prot1, col_ctr_Prot, col_exp_Prot)
f.row.stats <- function(a.df, a.col_ctr, a.col_exp) {
  a.df[, "ctr_count"] <- apply(!is.na(a.df[, a.col_ctr]), 1, sum)
  a.df[, "ctr_mean"]  <- rowMeans(a.df[, a.col_ctr], na.rm = TRUE)
  a.df[, "ctr_SD"]    <-
    rowSds(as.matrix(a.df[, a.col_ctr]), na.rm = TRUE)
  a.df[, "ctr_CV"]    <- a.df[, "ctr_SD"] / a.df[, "ctr_mean"] * 100
  
  a.df[, "exp_count"] <- apply(!is.na (a.df[, a.col_exp]), 1, sum)
  a.df[, "exp_mean"]  <- rowMeans(a.df[, a.col_exp], na.rm = TRUE)
  a.df[, "exp_SD"]    <-
    rowSds(as.matrix(a.df[, a.col_exp]), na.rm = TRUE)
  a.df[, "exp_CV"]    <- a.df[, "exp_SD"] / a.df[, "exp_mean"] * 100
  
  return(a.df)
}




# For LFQbenchmark experiments, perform a row of basic data sanitation
# to prepare for later following data filtering.
# This contains annotating entries with respective species 
# and retaining only single-species entries.
# e.g. Prot1 <- f.LFQbenchmark.basics(Prot0, col_ctr_Prot, col_exp_Prot)
f.LFQbenchmark.basics <- function(a.df,
                              a.col_ctr,
                              a.col_exp) {
  

  # Replace 0 with NA. "0" is sometimes used as mock value.
  # The mock value should not enter any calculation
  # and the event should be seen as missing value.(NA)
  a.df[a.col_ctr][a.df[a.col_ctr] == 0] <- NA
  a.df[a.col_exp][a.df[a.col_exp] == 0] <- NA
  
  
  
  # Calculate log2FC,
  # can be slightly different (5^-12 units) from limma,
  # but is easier to reproduce.
  a.df[, "log2FC"] <- log2(a.df[, "exp_mean"] / a.df[, "ctr_mean"])
  
  
  
  # In LFQbenchmark experiments, assign the respective species 
  # to species-specific entries and remove entries matching 
  # multiple species or no species. Should also work on 2-species.
  a.df[, "Species"] <- NA
  

  a.df[, "Species"][grepl('_HUMAN', a.df[, "Protein.Names"]) &
                      !grepl('_ECOLI', a.df[, "Protein.Names"]) &
                      !grepl('_YEAST', a.df[, "Protein.Names"]) &
                      !grepl('_CAEEL', a.df[, "Protein.Names"])] <- 'Human'
  
  a.df[, "Species"][grepl('_YEAST', a.df[, "Protein.Names"]) &
                      !grepl('_HUMAN', a.df[, "Protein.Names"]) &
                      !grepl('_ECOLI', a.df[, "Protein.Names"]) &
                      !grepl('_CAEEL', a.df[, "Protein.Names"])] <- 'Yeast'
  
  a.df[, "Species"][grepl('_ECOLI', a.df[, "Protein.Names"]) &
                      !grepl('_HUMAN', a.df[, "Protein.Names"]) &
                      !grepl('_YEAST', a.df[, "Protein.Names"]) &
                      !grepl('_CAEEL', a.df[, "Protein.Names"])] <- 'E.coli'
  
  a.df[, "Species"][grepl('_CAEEL', a.df[, "Protein.Names"]) &
                      !grepl('_HUMAN', a.df[, "Protein.Names"]) &
                      !grepl('_YEAST', a.df[, "Protein.Names"]) &
                      !grepl('_ECOLI', a.df[, "Protein.Names"])] <- 'C.elegans'
  

  a.df <- a.df[complete.cases(a.df[, "Species"]), ]
  
  
  
  
  
  # Add column listing expected log2FC row-wise
  # to simplify calculations of summary_stats around the log2FC.

    a.df[, "expected_log2FC"] <- NA
  

    a.df[a.df[, "Species"]== "Human", "expected_log2FC"]      <- expFC_human
    a.df[a.df[, "Species"]== "Yeast", "expected_log2FC"]      <- expFC_yeast
    a.df[a.df[, "Species"]== "E.coli", "expected_log2FC"]     <- expFC_ecoli
    a.df[a.df[, "Species"]== "C.elegans", "expected_log2FC"]  <- expFC_celegans
  
  
  return(a.df)
  
}




# Acquire adjusted p-values for differential expression analysis 
# for the protein group matrix with limma. Not executed on Precursor level.
# e.g. Prot3 <- f.limma.p_adj(Prot3, col_ctr_Prot, col_exp_Prot)
f.limma.p_adj <- function(a.df,
                          a.col_ctr,
                          a.col_exp) {
  
  # Limma to get BH adjusted p-values (p_adj) for differential expression.
  dcontra = matrix(c('exp', 'con'),
                   ncol = 2,
                   dimnames = list(NULL, c('CONDITION', 'CONDITION_REF'))) # defines in each rwo conditions that should be compared
  assay = as.matrix(a.df[, c(a.col_exp, a.col_ctr)]) # features x samples: Protein.Group x exp and con columns
  rownames(assay) = a.df[, 'Protein.Group']
  
  
 ## IF e.g. grps=c('exp','exp','exp','ctr','ctr','ctr')
  ## then for paired analysis use reps=c(1,2,3,1,2,3)
  ## Otherwise "reps" and "limma.rep" lines should be excluded from running.
  
    grps = c(rep('exp', length(a.col_exp)), rep('con', length(a.col_ctr))) # conditions of samples in assay
  
  
  # in case of paired analysis e.g. reps=c(1,2,3,1,2,3) 
  # reps = c(1:length(a.col_exp), 1:length(a.col_ctr))
  
  
  # Limma expects logarithmized input intensities
  assay = log2(assay)
  
  limma_data = ExpressionSet(assay
                             ,
                             phenoData = AnnotatedDataFrame(data.frame(
                               CONDITION = grps
                               , row.names = colnames(assay)
                             ))
                             ,
                             featureData = AnnotatedDataFrame(
                               data.frame(
                                 PROTEIN_IDS = a.df[, 'Protein.Ids']
                                 ,
                                 PROTEIN_NAMES = a.df[, 'Protein.Names']
                                 ,
                                 GENES = a.df[, 'Genes']
                                 ,
                                 row.names = rownames(assay)
                               )
                             ))
  comparisons <-
    c('exp - con') # Conditions that should be compared, in the trivial case it's exp vs con
  limma.cond = factor(grps)
  
  ## in case of paired analysis enable limma.rep
  #limma.rep <- factor(reps)
  contrast.matrix <- model.matrix(~  0 + limma.cond
                                 # + limma.rep
                                  )
  colnames(contrast.matrix) <-
    gsub("limma.cond", "", colnames(contrast.matrix))
  
  limma.object <- eBayes(contrasts.fit(
    lmFit(limma_data, design = contrast.matrix)
    ,
    makeContrasts(contrasts = comparisons, levels = contrast.matrix)
  )
  ,
  trend = TRUE,
  robust = TRUE)
  res_limma = limma::topTable(
    limma.object,
    coef = comparisons[1],
    number = Inf,
    sort.by = 'none',
    adjust.method = "BH",
    confint = TRUE
  )
  
  
  # Add p_adj (limma result) to the main protein group data frame.
  a.df[, "p_adj"] <- res_limma[, 'adj.P.Val']
  a.df[, "p_adj"] <- as.numeric(format(a.df[, "p_adj"],
                                       scientific = FALSE,
                                       justified = "none"))
  return(a.df)
  
  
  # # For use outside of a function,
  # # clean up environment from unused limma elements.
  # rm(assay)
  # rm(contrast.matrix)
  # rm(dcontra)
  # rm(limma_data)
  # rm(limma.object)
  # rm(comparisons)
  # rm(grps)
  # rm(limma.cond)
  # rm(limma.rep)
  # rm(reps)
  # rm(res_limma)
}




# "Interpret" the p_adj from limma differential expression analysis,
# this generates the "confusion matrix" summary stats.
# !!! Hard-coded to the orientation of the example in the introduction.
# Please see introduction.
# e.g. Prot3 <- f.LFQbenchmark.limma.interpretation(Prot3, alpha_limma, limit_FC)
f.LFQbenchmark.limma.interpretation <-
  function(a.df, a.alpha_limma, a.limit_FC) {
    
    a.df[(a.df[, "Species"] == "Human"
          & (a.df[, "p_adj"] >= a.alpha_limma
             | abs(a.df[, "log2FC"]) <= a.limit_FC)), "DE_result"] <-
      "true negative"
    
    a.df[(a.df[, "Species"] == "Human"
          & a.df[, "p_adj"] < a.alpha_limma
          & abs(a.df[, "log2FC"]) > a.limit_FC), "DE_result"] <-
      "false positive"
    
    
    
    a.df[(a.df[, "Species"] == "Yeast"
          & a.df[, "log2FC"] > +a.limit_FC
          & a.df[, "p_adj"] < a.alpha_limma), "DE_result"] <-
      "true positive"
    
    a.df[(a.df[, "Species"] == "Yeast"
          & a.df[, "log2FC"] < -a.limit_FC
          & a.df[, "p_adj"] < a.alpha_limma), "DE_result"] <-
      "false positive"
    
    a.df[(a.df[, "Species"] == "Yeast"
          & (a.df[, "p_adj"] >= a.alpha_limma
             | abs(a.df[, "log2FC"]) <= a.limit_FC)), "DE_result"] <-
      "false negative"
    
    
    
    a.df[(a.df[, "Species"] == "E.coli"
          & a.df[, "log2FC"] < -a.limit_FC
          & a.df[, "p_adj"] < a.alpha_limma), "DE_result"] <-
      "true positive"
    
    a.df[(a.df[, "Species"] == "E.coli"
          & a.df[, "log2FC"] > +a.limit_FC
          & a.df[, "p_adj"] < a.alpha_limma), "DE_result"] <-
      "false positive"
    
    a.df[(a.df[, "Species"] == "E.coli"
          & (a.df[, "p_adj"] >= a.alpha_limma
             | abs(a.df[, "log2FC"]) <= a.limit_FC)), "DE_result"] <-
      "false negative"
    
    
    
    a.df[(a.df[, "Species"] == "C.elegans"
          & a.df[, "log2FC"] < -a.limit_FC
          & a.df[, "p_adj"] < a.alpha_limma), "DE_result"] <-
      "true positive"
    
    a.df[(a.df[, "Species"] == "C.elegans"
          & a.df[, "log2FC"] > +a.limit_FC
          & a.df[, "p_adj"] < a.alpha_limma), "DE_result"] <-
      "false positive"
    
    a.df[(a.df[, "Species"] == "C.elegans"
          & (a.df[, "p_adj"] >= a.alpha_limma
             | abs(a.df[, "log2FC"]) <= a.limit_FC)), "DE_result"] <-
      "false negative"
    
    
    return(a.df)
    
  }




# Tailing and asymmetry factor for up- and down-regulated proteins (according to species mixtures)
# Different calculations to acquire tailing towards the outside away from log2FC of 0 (or vice versa),
# in both directions.
# E.g. f.asym.down-regulated(Prot3, "E.coli")
f.asym.downregulated <- function(a.df,
                                 a.species) {
  
  # if statement required to prevent errors when handling 2-species benchmarks.
  if(length(a.df[a.df[, "Species"] == a.species, "log2FC"]) > 2) {
    
    # 1-dimensional data like log2 fold-changes
    # to data frames containing x and y coordinates of density plot function
    data <- a.df[a.df[, "Species"] == a.species, "log2FC"]
    density <- density(data)
    density <- as.data.frame(cbind(density$x, density$y))
    colnames(density) <- c("x", "y")
    
    # ymax and respective x coordinate called C
    ymax <- density$y[which.max(density$y)]
    C <- density$x[which.max(density$y)]
    
    # subset density coordinates by lefthand and righthand of modal
    # to differentiate multiple x coordinates matching one y coordinate.
    density_left <- density[density$x < C,] 
    density_right <- density[density$x > C,] 
    
    # A and B are the index, x, and y coordinates at a given height relative to ymax.
    A1 <- density_left[which.min(abs(density_left$y - ymax*0.10)),]
    B1 <- density_right[which.min(abs(density_right$y - ymax*0.10)),]
    
    A2 <- density_left[which.min(abs(density_left$y - ymax*0.05)),]
    B2 <- density_right[which.min(abs(density_right$y - ymax*0.05)),]
    
    # Factors calculated directional towards the outside away from log2FC of 0.
    Asymmetry_Factor <- (A1$x - C) / (C - B1$x)
    Tailing_Factor <- (A2$x - B2$x) / (2* (C - B2$x))
    
    return(cbind(Tailing_Factor, Asymmetry_Factor))
    
  } else {
    
    return(cbind(NA, NA))}
  
}

# Asymmetry and tailing calculated for up-regulated log2 fold-changes.
# # E.g. f.asym.up-regulated(Prot3, "Yeast")
f.asym.upregulated <- function(a.df,
                               a.species) {
  
  # if statement required to prevent errors when handling 2-species benchmarks.
  if(length(a.df[a.df[, "Species"] == a.species, "log2FC"]) > 2) {
    
    # 1-dimensional data like log2 fold-changes
    # to data frames containing x and y coordinates of density plot function
    data <- a.df[a.df[, "Species"] == a.species, "log2FC"]
    density <- density(data)
    density <- as.data.frame(cbind(density$x, density$y))
    colnames(density) <- c("x", "y")
    
    # ymax and respective x coordinate called C
    ymax <- density$y[which.max(density$y)]
    C <- density$x[which.max(density$y)]
    
    # subset density coordinates by leftside and rightside of modal
    # to differentiate x coordinates matchign one y coordiante
    density_left <- density[density$x < C,] 
    density_right <- density[density$x > C,] 
    
    # A and B are the index, x, and y coordinates at a given height relative to ymax.
    A1 <- density_left[which.min(abs(density_left$y - ymax*0.10)),]
    B1 <- density_right[which.min(abs(density_right$y - ymax*0.10)),]
    
    A2 <- density_left[which.min(abs(density_left$y - ymax*0.05)),]
    B2 <- density_right[which.min(abs(density_right$y - ymax*0.05)),]
    
    # Factors calculated directional towards the outside away from log2FC of 0.
    Asymmetry_Factor <- (B1$x - C) / (C - A1$x)
    Tailing_Factor <- (B2$x - A2$x) / (2* (C - A2$x))
    
    return(cbind(Tailing_Factor, Asymmetry_Factor))
    
  } else {
    
    return(cbind(NA, NA))}
  
}




# For processed LFQbenchmark data, collect main summary statistics
# into a new data frame for export and plot subtitles.
# They describe sensitivity, precision, accuracy, reliability, etc..
# Next to plots they are the core result of this script.
# Mainly use true positives (TP), deFDR, and Asymmetry_Factor to evaluate the performance.
# See introduction for interpretation.
f.LFQbenchmark.summary.stats <-
  function(a.df01, a.df02, a.df03, a.df04) {
    summary_stats <- data.frame(matrix(ncol = 0, nrow = 1))
    
    
    # R script filter variables to avoid confusion
    # when collecting stats of multiple benchmarks in excel.
    summary_stats[, "Variables"] <-
      paste0(script_version,",",
             rep_min, "of", rep_max,",",
             limit_CV,",",
             limit_FC,",",
             alpha_limma)
    
    
    # deFDR (differential expression FDR) (part of confusion matrix)
    summary_stats[, "deFDR"] <-
      round(digits = 2,
            100 * length(which(a.df02[, "DE_result"] == "false positive")) /
              (length(which(
                a.df02[, "DE_result"] == "false positive"
              ))
              + length(which(
                a.df02[, "DE_result"] == "true positive"
              ))))
    
    
    # true positive count
    summary_stats[, "TP"] <-
      length(which(a.df02[, "DE_result"] == "true positive"))
    
    
    # Number of confidently detected and quantified protein groups,
    # after filtering for missingness and missingness + CV.
    summary_stats[, "Prot_ID"] <- nrow(a.df01)
    summary_stats[, "Prot_Quant"] <- nrow(a.df02)
    
    
    # Asymmetry_Factor to evaluate ratio compression or extension.
    summary_stats[, "Prot_Asymmetry_E.coli"] <- 
      round( digits =2,
             as.numeric(subset(asymmetry, 
                               Species == "E.coli" &
                                 Group == "protein group",
                               select = "Asymmetry_Factor")))
    
    summary_stats[, "Prot_Asymmetry_Yeast"] <- 
      round( digits =2,
             as.numeric(subset(asymmetry, 
                               Species == "Yeast" &
                                 Group == "protein group",
                               select = "Asymmetry_Factor")))
    
    
    # Coefficient of variation to describe precision.
    # For ca. 3 replicates aim at average and median at or below 5%
    # after filtering for <20%.
    summary_stats[, "Prot_CV_Mean"]   <-
      round(digits = 2,  mean(c(a.df02[, "ctr_CV"], a.df02[, "exp_CV"])))
    summary_stats[, "Prot_CV_Median"] <-
      round(digits = 2, median(c(a.df02[, "ctr_CV"], a.df02[, "exp_CV"])))
    
    
    # Accuracy is average distance of measured to expected log2FC.
    summary_stats[, "Prot_Accuracy"] <-
      round(digits = 2,
            mean(abs(a.df02[, "log2FC"] - a.df02[, "expected_log2FC"])))
    
    
    # Trueness is distance between measurement medians 
    # and respective expected fold-changes.
    # (cumulative as shifts are important to spot)
    # Several conditions can lead to a high value
    # e.g. ratio compression or erroneous normalization or even both together.
    summary_stats[, "Prot_Trueness"] <-
      round(
        digits = 2,
        sum(na.rm = TRUE,
            c(abs(median(a.df02[a.df02[, "Species"] == "Human"  , "log2FC"]) - expFC_human),
              abs(median(a.df02[a.df02[, "Species"] == "Yeast"  , "log2FC"]) - expFC_yeast),
              abs(median(a.df02[a.df02[, "Species"] == "E.coli" , "log2FC"]) - expFC_ecoli),
              abs(median(a.df02[a.df02[, "Species"] == "C.elegans"  , "log2FC"]) - expFC_celegans)
              )))
    
    
    # Dispersion is the average distance of individual measurements around the respective median.
    summary_stats[, "Prot_Dispersion"] <-
      round(digits = 2,
            mean(na.rm = TRUE,
                 c(
                   abs(a.df02[a.df02[, "Species"] == "Human"   , "log2FC"] - median(a.df02[a.df02[, "Species"] == "Human"  , "log2FC"])),
                   abs(a.df02[a.df02[, "Species"] == "Yeast"   , "log2FC"] - median(a.df02[a.df02[, "Species"] == "Yeast"  , "log2FC"])),
                   abs(a.df02[a.df02[, "Species"] == "E.coli"  , "log2FC"] - median(a.df02[a.df02[, "Species"] == "E.coli" , "log2FC"])),
                   abs(a.df02[a.df02[, "Species"] == "C.elegans"  , "log2FC"] - median(a.df02[a.df02[, "Species"] == "C.elegans" , "log2FC"]))
                 )))
    
    
    # Other than TP, continued confusion matrix of 
    # differential expression result interpretation.
    summary_stats[, "FP"] <-
      length(which(a.df02[, "DE_result"] == "false positive"))
    summary_stats[, "TN"] <-
      length(which(a.df02[, "DE_result"] == "true negative"))
    summary_stats[, "FN"] <-
      length(which(a.df02[, "DE_result"] == "false negative"))
    
    
    summary_stats[, "Sensitivity"] <-
      round(digits = 2, 100 * length(which(a.df02[, "DE_result"] == "true positive")) /
              (length(which(
                a.df02[, "DE_result"] == "true positive"
              )) + length(which(
                a.df02[, "DE_result"] == "false negative"
              ))))
    
    summary_stats[, "Specificity"] <-
      round(digits = 2, 100 * length(which(a.df02[, "DE_result"] == "true negative")) /
              (length(which(
                a.df02[, "DE_result"] == "true negative"
              )) + length(which(
                a.df02[, "DE_result"] == "false positive"
              ))))
    
    
    
    
    
    
    # Stats partially repeated on precursor level as for protein groups above.
    summary_stats[, "Prec_ID"] <- nrow(a.df03)
    summary_stats[, "Prec_Quant"] <- nrow(a.df04)
    
    
    # Asymmetry_Factor to evaluate ratio compression or extension.
    summary_stats[, "Prec_Asymmetry_E.coli"] <- 
      round( digits =2,
             as.numeric(subset(asymmetry, 
                               Species == "E.coli" &
                                 Group == "precursor",
                               select = "Asymmetry_Factor")))
    
    summary_stats[, "Prec_Asymmetry_Yeast"] <- 
      round( digits =2,
             as.numeric(subset(asymmetry, 
                               Species == "Yeast" &
                                 Group == "precursor",
                               select = "Asymmetry_Factor")))
    
    
    
    # Precursor-level precision by coefficient of variation.
    summary_stats[, "Prec_CV_Mean"]   <-  round(digits = 2,
                                                mean(c(a.df04[, "ctr_CV"], a.df04[, "exp_CV"])))
    summary_stats[, "Prec_CV_Median"] <-  round(digits = 2,
                                                median(c(a.df04[, "ctr_CV"], a.df04[, "exp_CV"])))
    
    
    # Accuracy is average distance of measured to expected log2FC.
    summary_stats[, "Prec_Accuracy"] <-
      round(digits = 2,
            mean(abs(a.df04[, "log2FC"] - a.df04[, "expected_log2FC"])))
    
    
    # Trueness is distance between measurement medians 
    # and respective expected fold-changes (cumulative as shifts are important to spot)
    summary_stats[, "Prec_Trueness"] <-
      round(
        digits = 2,
        sum(na.rm = TRUE,
            c(abs(median(a.df04[a.df04[, "Species"] == "Human"  , "log2FC"]) - expFC_human),
              abs(median(a.df04[a.df04[, "Species"] == "Yeast"  , "log2FC"]) - expFC_yeast),
              abs(median(a.df04[a.df04[, "Species"] == "E.coli"  , "log2FC"]) - expFC_ecoli),
              abs(median(a.df04[a.df04[, "Species"] == "C.elegans"  , "log2FC"]) - expFC_celegans)
              )))
    
    
    # Dispersion is the average distance of individual measurements around the respective median.
    summary_stats[, "Prec_Dispersion"] <-
      round(digits = 2,
            mean(na.rm = TRUE,
                 c(
                   abs(a.df04[a.df04[, "Species"] == "Human"  , "log2FC"] - median(a.df04[a.df04[, "Species"] == "Human"  , "log2FC"])),
                   abs(a.df04[a.df04[, "Species"] == "Yeast"  , "log2FC"] - median(a.df04[a.df04[, "Species"] == "Yeast"  , "log2FC"])),
                   abs(a.df04[a.df04[, "Species"] == "E.coli"  , "log2FC"] - median(a.df04[a.df04[, "Species"] == "E.coli"  , "log2FC"])),
                   abs(a.df04[a.df04[, "Species"] == "C.elegans"  , "log2FC"] - median(a.df04[a.df04[, "Species"] == "C.elegans"  , "log2FC"]))
                 )))
    
    
    
    # Folder that contained the analysed pg and pc matrix
    summary_stats[, "folder_input"] <- folder_input
    
    
    
    
    return(summary_stats)
    
  }





# Import, Processing, Filtering ------------------------------------------------------------------


# Import protein group and precursor matrix.
Prot0 <- f.import.matrix(folder_input, "pg_matrix", "\t")
Prec0 <- f.import.matrix(folder_input, "pr_matrix", "\t")



# Perform basic sanitation and calculations as preparation, 
# but no quantification-related filtering.
Prot1 <- Prot0
Prot1 %<>%
  f.strip.sample.column.name(c(col_ctr_Prot, col_exp_Prot)) %>%
  f.remove.MQ.cont() %>%
  f.row.stats(col_ctr_Prot, col_exp_Prot) %>%
  f.LFQbenchmark.basics(col_ctr_Prot, col_exp_Prot)


# Retain entries matching defined data missingness or completeness.
# Passing entries count as "IDs".
Prot2 <- subset(Prot1, ctr_count >= rep_min & exp_count >= rep_min)


# Continue filtering by precision (CV)
# Passing entries count as "Quant." (quantifications)
Prot3 <- subset(Prot2, ctr_CV <  limit_CV & exp_CV <  limit_CV)


# Differential expression analysis and match of results
# against expectations from sample mixtures.
Prot3 %<>%
  f.limma.p_adj(col_ctr_Prot, col_exp_Prot) %>%
  f.LFQbenchmark.limma.interpretation(alpha_limma, limit_FC)






# Repeat data processing except for differential expression
# analysis on the precursor level.

# Perform basic sanitation and calculations as preparation, 
# but no quantification-related filtering.
Prec1 <- Prec0
Prec1 %<>%
  f.strip.sample.column.name(c(col_ctr_Prec, col_exp_Prec)) %>%
  f.remove.MQ.cont() %>%
  f.row.stats(col_ctr_Prec, col_exp_Prec) %>%
  f.LFQbenchmark.basics(col_ctr_Prec, col_exp_Prec)


# Retain entries matching defined data missingness or completeness.
# Passing entries count as "IDs".
Prec2 <- subset(Prec1, ctr_count >= rep_min & exp_count >= rep_min)


# Continue filtering by precision (CV)
# Passing entries count as "Quant." (quantifications)
Prec3 <- subset(Prec2, ctr_CV <  limit_CV & exp_CV <  limit_CV)





# Summary Stats -----------------------------------

# At first produce a separate data frame and output
# for summary stats describing tendencies of fold-change over-or underestimation,
# then the main data frame and export of overall main summary_stat.



# Summary stats describing potential asymmetry in density plots
# to investigate fold-change over-or underestimation.
# Quantiles will be attached below, but the asymmetry factor remains
# the easiest to interpret and the most important of these summary stats.
# Asymmetry Factor of 0.5 or below indicates a strong fold-change underestimation, 
# a value above 2 indicates a strong fold-change overestimation.

asymmetry <- data.frame(matrix(ncol = 0, nrow = 8))
asymmetry[, "Group"] <- c(rep("precursor", 4), 
                          rep("protein group", 4))
asymmetry[, "Species"] <- c("E.coli", "Human", "Yeast", "C.elegans")



asymmetry[, c("Asymmetry_Factor")] <- c(f.asym.downregulated(Prec3, "E.coli")[,2],
                                        f.asym.downregulated(Prec3, "Human")[,2],
                                        f.asym.upregulated(Prec3, "Yeast")[,2],
                                        f.asym.upregulated(Prec3, "C.elegans")[,2],
                                        f.asym.downregulated(Prot3, "E.coli")[,2],
                                        f.asym.downregulated(Prot3, "Human")[,2],
                                        f.asym.upregulated(Prot3, "Yeast")[,2],
                                        f.asym.upregulated(Prot3, "C.elegans")[,2])

asymmetry[, c("Tailing_Factor")] <- c(f.asym.downregulated(Prec3, "E.coli")[,1],
                                      f.asym.downregulated(Prec3, "Human")[,1],
                                      f.asym.upregulated(Prec3, "Yeast")[,1],
                                      f.asym.upregulated(Prec3, "C.elegans")[,1],
                                      f.asym.downregulated(Prot3, "E.coli")[,1],
                                      f.asym.downregulated(Prot3, "Human")[,1],
                                      f.asym.upregulated(Prot3, "Yeast")[,1],
                                      f.asym.upregulated(Prot3, "C.elegans")[,1])



asymmetry[, "Skewness"] <- c(
  skewness(Prec3[which(Prec3[, "Species"] == "E.coli"), "log2FC"]),
  skewness(Prec3[which(Prec3[, "Species"] == "Human"), "log2FC"]),
  skewness(Prec3[which(Prec3[, "Species"] == "Yeast"), "log2FC"]),
  skewness(Prec3[which(Prec3[, "Species"] == "C.elegans"), "log2FC"]),
  skewness(Prot3[which(Prot3[, "Species"] == "E.coli"), "log2FC"]),
  skewness(Prot3[which(Prot3[, "Species"] == "Human"), "log2FC"]),
  skewness(Prot3[which(Prot3[, "Species"] == "Yeast"), "log2FC"]),
  skewness(Prot3[which(Prot3[, "Species"] == "C.elegans"), "log2FC"]))

asymmetry[, "Kurtosis"] <- c(
  kurtosis(Prec3[which(Prec3[, "Species"] == "E.coli"), "log2FC"]),
  kurtosis(Prec3[which(Prec3[, "Species"] == "Human"), "log2FC"]),
  kurtosis(Prec3[which(Prec3[, "Species"] == "Yeast"), "log2FC"]),
  kurtosis(Prec3[which(Prec3[, "Species"] == "C.elegans"), "log2FC"]),
  kurtosis(Prot3[which(Prot3[, "Species"] == "E.coli"), "log2FC"]),
  kurtosis(Prot3[which(Prot3[, "Species"] == "Human"), "log2FC"]),
  kurtosis(Prot3[which(Prot3[, "Species"] == "Yeast"), "log2FC"]),
  kurtosis(Prot3[which(Prot3[, "Species"] == "C.elegans"), "log2FC"]))


# Protein group-level quantiles.
quantiles_Prot <- Prot3 %>%
  group_by(Species) %>%
  summarise(Q1 = quantile(log2FC, probs = 0.25),
            Q2 = quantile(log2FC, probs = 0.50),
            Q3 = quantile(log2FC, probs = 0.75))

quantiles_Prot[,"Group"] <- "protein group"
quantiles_Prot[,"|Q2-Q1|"] <- abs(quantiles_Prot[,"Q2"] - quantiles_Prot[,"Q1"])
quantiles_Prot[,"|Q2-Q3|"] <- abs(quantiles_Prot[,"Q2"] - quantiles_Prot[,"Q3"])
quantiles_Prot[,"|Q2-Q1| / |Q2-Q3|"] <- quantiles_Prot[,"|Q2-Q1|"] / quantiles_Prot[,"|Q2-Q3|"] 
quantiles_Prot[,"|Q2-Q3| / |Q2-Q1|"] <- quantiles_Prot[,"|Q2-Q3|"] / quantiles_Prot[,"|Q2-Q1|"]


# Precursor level quantiles.
quantiles_Prec <- Prec3 %>%
  group_by(Species) %>%
  summarise(Q1 = quantile(log2FC, probs = 0.25),
            Q2 = quantile(log2FC, probs = 0.50),
            Q3 = quantile(log2FC, probs = 0.75))

quantiles_Prec[,"Group"] <- "precursor"
quantiles_Prec[,"|Q2-Q1|"] <- abs(quantiles_Prec[,"Q2"] - quantiles_Prec[,"Q1"])
quantiles_Prec[,"|Q2-Q3|"] <- abs(quantiles_Prec[,"Q2"] - quantiles_Prec[,"Q3"])
quantiles_Prec[,"|Q2-Q1| / |Q2-Q3|"] <- quantiles_Prec[,"|Q2-Q1|"] / quantiles_Prec[,"|Q2-Q3|"] 
quantiles_Prec[,"|Q2-Q3| / |Q2-Q1|"] <- quantiles_Prec[,"|Q2-Q3|"] / quantiles_Prec[,"|Q2-Q1|"]



# Bind precursor and protein group level quantiles.
quantiles <- rbind( quantiles_Prec, quantiles_Prot)
rm(quantiles_Prot,quantiles_Prec)


# Bind asymmetry stats and quantiles to form full overview on asymmetry stats.
asymmetry <- merge(asymmetry, quantiles,
           by = c("Group", "Species"))




# Export combined asymmetry-related summary stats.
write.csv(asymmetry,
          row.names = FALSE,
          paste0(folder_output, "/", "asymmetry_stats_",
                 basename(dirname(folder_output)),"_",
                 script_version,"_",
                 rep_min,"of", rep_max,
                 "_CV", limit_CV,
                 "_FC",limit_FC,
                 "_a",alpha_limma, ".csv"))



# Main summary stats.




# Create data frame of LFQbenchmark main summary stats (main result).
summary_stats <- f.LFQbenchmark.summary.stats(Prot2, Prot3, Prec2, Prec3)

# Export summary stats.
write.csv(summary_stats,
          row.names = FALSE,
          paste0(folder_output, "/", "summary_stats_",
                 basename(dirname(folder_output)),"_",
                 script_version,"_",
                 rep_min,"of", rep_max,
                 "_CV", limit_CV,
                 "_FC",limit_FC,
                 "_a",alpha_limma, ".csv"))




# Plot preparations -------------------------------------------------------

# Title, just pretyped.
title_Prot <- NA
title_Prec <- title_Prot


# Create subtitles with selection of summary stats.
subtitle_Prot <- paste0(
  summary_stats[, "Prot_Quant"],
  " Quant. (",
  round(nrow(Prot3) / nrow(Prot2) * 100, digits = 0),
  "%), ",
  summary_stats[, "Prot_ID"],
  " IDs",
  " \n Accu ",
  summary_stats[, "Prot_Accuracy"],
  ", ",
  "True ",
  summary_stats[, "Prot_Trueness"],
  ", ",
  "Disp ",
  summary_stats[, "Prot_Dispersion"],
  "\n",
  summary_stats[, "TP"],
  " true positives, ",
  summary_stats[, "deFDR"],
  "% deFDR")



  subtitle_Prec <- paste0(
    summary_stats[, "Prec_Quant"],
    " Quant. (",
    round(nrow(Prec3) / nrow(Prec2) * 100, digits = 0),
    "%), ",
    summary_stats[, "Prec_ID"],
    " IDs",
    " \n Accu ",
    summary_stats[, "Prec_Accuracy"],
    ", ",
    "True ",
    summary_stats[, "Prec_Trueness"],
    ", ",
    "Disp ",
    summary_stats[, "Prec_Dispersion"],
    "\n ")
  

# Different captions for precursor and protein group plots.
caption_Prot <- paste0(rep_min, "of", rep_max,
                       ", CV < ", limit_CV, "%, ",
                       "|log2FC| > ", limit_FC)

caption_Prec <- paste0(rep_min, "of", rep_max,
                       ", CV < ", limit_CV, "%")


# Plotting order of species.
Prot3$Species <-
  factor(Prot3$Species, levels = c("Human", "Yeast", "E.coli", "C.elegans"))
Prec3$Species <-
  factor(Prec3$Species, levels = c("Human", "Yeast", "E.coli", "C.elegans"))






# Function to export plots to folder_output.
f.export.plot <- function(a.plotname,
                          a.filename,
                          a.width,
                          a.height,
                          a.res) {
  png(
    file = paste0(
      folder_output,
      "/",
      a.filename,
      "_",
      basename(dirname(folder_output)),
      ".png"
    ),
    units = "cm",
    width = a.width,
    height = a.height,
    res = a.res
  )
  print(a.plotname)
  dev.off()
}






# Plot Replicate PCA ---------------------------------------------------------------------

f.p.pca <- function(a.df,
                    a.col_ctr,
                    a.col_exp,
                    a.subtitle,
                    a.caption) {

  
# Replicate-level pca calculations
pca <- prcomp(t(na.omit(a.df[,c(a.col_ctr, a.col_exp)])), scale = TRUE)

pca_df <- data.frame(Replicate = rownames(pca$x),
                       x = pca$x[,1],
                       y = pca$x[,2])

pca_var <- round(pca$sdev^2/sum(pca$sdev^2)*100,1)


# Replicate pca plot
ggplot(pca_df,
       aes(
         x = x,
         y = y,
         label = Replicate)) +
  
  coord_cartesian(
    xlim = c(2.5*min(pca_df[,"x"]),2.5*max(pca_df[,"x"])),
    ylim = c(2.5*min(pca_df[,"y"]),2.5*max(pca_df[,"y"])))+
  
  geom_point(size = 2.5, alpha = 0.5) +
  
  geom_text(size = 2.0, 
            vjust= 1.9)+
  
  labs(subtitle = a.subtitle,
       caption = a.caption)+
  
       xlab(paste0("PC1 - ", pca_var[1],"%")) +
       ylab(paste0("PC2 - ", pca_var[2],"%")) +
  

  theme(
    panel.background  = element_blank(),
    plot.background = element_rect(fill = "white", color = "white"),
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      color = 'black',
      size = 12,
      angle = 0,
      margin = margin(0, 0, 0.0, 0)
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5,
      face = "plain",
      color = 'black',
      size = 12,
      angle = 0,
      margin = margin(0, 0, 1, 0)
    ),
    
    plot.caption = element_text(
      hjust = -0.1,
      face = "italic",
      color = 'black',
      size = 7.5,
      angle = 0,
      margin = margin(5, 0,-4, 0)
    ),
    
    axis.text.x = element_text(
      face = "plain",
      color = 'black',
      size = 10,
      angle = 0,
      margin = margin(2, 0, 0, 0)
    ),
    
    axis.text.y = element_text(
      face = "plain",
      color = 'black',
      size = 10,
      angle = 0,
      margin = margin(0, 2, 0, 0)
    ),
    
    axis.title.y = element_text(
      face = "plain",
      color = 'black',
      size = 11,
      angle = 90,
      margin = margin(0, 4, 0, 0)
    ),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.line = element_line(
      size = 0.3,
      linetype = "solid",
      color = "black"
    ),
    
    axis.ticks = element_line(
      size = 0.3,
      linetype = "solid",
      color = "black"
    ),
    
    axis.title.x = element_text(
      hjust = 0.5,
      face = "plain",
      color = 'black',
      size = 11,
      angle = 0,
      margin = margin(4, 0, 0, 0)
    ),
    
    # legend.title = element_blank(),
    # legend.background = element_rect(fill = "transparent",
    #                                  color = NA, size = 0.1),
    #
    legend.key = element_rect(fill = "transparent", color = NA),
    # legend.key.height = unit(0, "cm"),
    # legend.margin = margin(3, 0, 0, 0),
    # legend.text.align = 0,
    # legend.box.spacing = unit(0.0, "cm"),
    # legend.box.margin = margin(0, 0, 0, -0.55, "cm"),
    # legend.spacing.x = unit(0.0, "cm"),
    # legend.position = "top",
    # legend.text = element_text(face = "plain", color = 'black',
    #                            size = 9, angle = 0,
    #                            margin = margin(0, 0, 0, 0))
  ) +
  
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

}

# Generate plots.
plot_Prot_pca <-  f.p.pca(Prot3,
                          col_ctr_Prot,
                          col_exp_Prot,
                          subtitle_Prot,
                          caption_Prot)

plot_Prec_pca <-  f.p.pca(Prec3,
                          col_ctr_Prec,
                          col_exp_Prec,
                          subtitle_Prec,
                          caption_Prec)

# Export plots.
f.export.plot(plot_Prot_pca, "Prot_pca", plot_width, plot_height, plot_res)
f.export.plot(plot_Prec_pca, "Prec_pca", plot_width, plot_height, plot_res)





# Plot CV violin ----------------------------------------------------------

f.p.CV <- function(a.df,
                   a.subtitle,
                   a.caption) {
  ggplot((
    pivot_longer(
      a.df[, c("ctr_CV", "exp_CV", "Species")],
      cols = c("ctr_CV", "exp_CV"),
      names_to = "TechRep",
      values_to = "CV"
    )
  ),
  aes(x = Species,
      y = CV,
      fill = Species)) +
    
    geom_violin(
      trim = TRUE,
      scale = "area",
      col = "black",
      key_glyph = "blank",
      alpha = 0.6,
      draw_quantiles = c(0.5)
    ) +
    
    stat_summary(
      fun = mean,
      geom = "point",
      size = 1.5,
      color = "black",
      show.legend = FALSE
    ) +
    
    labs(subtitle = a.subtitle,
         caption = a.caption) +
    
    scale_fill_manual(values = c(color_human, color_yeast, color_ecoli, color_celegans)) +
    
    coord_cartesian(ylim = c(0, 30)) +
    scale_y_continuous(
      name = "CV in %",
      breaks = c(0, 5, 10, 15, 20, 25, 30),
      expand = expansion(mult = c(.00, .00))
    ) +
    
    theme_gray(base_size = 11) +
    theme(
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      plot.background = element_rect(fill = "white",
                                     color = "white"),
      
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        color = 'black',
        size = 12,
        angle = 0,
        margin = margin(0, 0, 0.0, 0)
      ),
      
      plot.subtitle = element_text(
        hjust = 0.5,
        face = "plain",
        color = 'black',
        size = 12,
        angle = 0,
        margin = margin(0, 0, 1, 0)
      ),
      
      plot.caption = element_text(
        hjust = -0.1,
        face = "italic",
        color = 'black',
        size = 7.5,
        angle = 0,
        margin = margin(5, 0,-4, 0)
      ),
      
      axis.text.x = element_text(
        face = "plain",
        color = 'black',
        size = 10,
        angle = 0,
        margin = margin(2, 0, 0, 0)
      ),
      
      axis.text.y = element_text(
        face = "plain",
        color = 'black',
        size = 10,
        angle = 0,
        margin = margin(0, 2, 0, 0)
      ),
      
      axis.title.x = element_text(
        hjust = 0.5,
        face = "plain",
        color = "white",
        size = 11,
        angle = 0,
        margin = margin(4, 0, 0, 0)
      ),
      
      axis.title.y = element_text(
        face = "plain",
        color = 'black',
        size = 11,
        angle = 90,
        margin = margin(0, 4, 0, 0)
      ),
      
      axis.line = element_line(
        size = 0.3,
        linetype = "solid",
        color = "black"
      ),
      
      axis.ticks = element_line(
        size = 0.3,
        linetype = "solid",
        color = "black"
      ),
      
      
      legend.background = element_rect(
        fill = "transparent",
        color = NA,
        size = 0.1
      ),
      legend.title = element_blank(),
      legend.key = element_rect(fill = "transparent", color = NA),
      legend.key.height = unit(0, "cm"),
      legend.margin = margin(2, 0, 0, 0),
      legend.text.align = 0,
      legend.box.spacing = unit(0.0, "cm"),
      legend.box.margin = margin(0, 0, 0,-0.55, "cm"),
      legend.spacing.x = unit(0.0, "cm"),
      legend.position = "top",
      legend.text = element_text(
        face = "plain",
        color = "white",
        size = 9,
        angle = 0,
        margin = margin(0, 0, 0, 0)
      )
    )
}


# Generate plots
plot_Prot_CV <- f.p.CV(Prot3, subtitle_Prot, caption_Prot)
plot_Prec_CV <- f.p.CV(Prec3, subtitle_Prec, caption_Prec)


# Export plots
f.export.plot(plot_Prot_CV, "Prot_CV", plot_width, plot_height, plot_res)
f.export.plot(plot_Prec_CV, "Prec_CV", plot_width, plot_height, plot_res)






# Plot density ------------------------------------------------------------

f.p.density <- function(a.df,
                        a.subtitle,
                        a.caption)   {
  ggplot(a.df,
         aes(x = log2FC,
             fill = factor(Species))) +
    
    # Values outside FC_min and FC_max are ignored,
    # to avoid coord_cartesian bug resulting in loss of plot color.
    scale_x_continuous(
      name = "log2(A:B)",
      breaks = seq(from = round(digits = 0, FC_min), to = round(digits = 0, FC_max), by = 1 ),
      expand = expansion(mult = c(.02, .00)),
      limits = c(FC_min, FC_max)
    ) +
    
    scale_y_continuous(
      name = "density",
      breaks = seq(from = 0, to = 7, by = 1 ),
      expand = expansion(mult = c(.00, .02)),
      limits = c(0, 7)
    ) +
    
    
    geom_vline(
      xintercept = median(a.df[a.df[, "Species"] == "Human", "log2FC"]),
      linetype = "dashed",
      color = color_human,
      alpha = 0.9
    ) +
    
    geom_vline(
      xintercept = median(a.df[a.df[, "Species"] == "Yeast", "log2FC"]),
      linetype = "dashed",
      color = color_yeast,
      alpha = 0.9
    ) +
    
    geom_vline(
      xintercept = median(a.df[a.df[, "Species"] == "E.coli", "log2FC"]),
      linetype = "dashed",
      color = color_ecoli,
      alpha = 0.9
    ) +
    
    geom_vline(
      xintercept = median(a.df[a.df[, "Species"] == "C.elegans", "log2FC"]),
      linetype = "dashed",
      color = color_celegans,
      alpha = 0.9
    ) +
    

    geom_density(
      data = a.df[a.df[, "Species"] == "Human", ],
      size = 0.6,
      alpha = 0.6,
      # key_glyph =  draw_key_rect,
      fill = color_human,
    ) +
    
    geom_density(
      data = a.df[a.df[, "Species"] == "Yeast", ],
      size = 0.6,
      alpha = 0.6,
      # key_glyph =  draw_key_rect,
      fill = color_yeast
      ) +
    
    geom_density(
      data = a.df[a.df[, "Species"] == "E.coli", ],
      size = 0.6,
      alpha = 0.6,
      # key_glyph =  draw_key_rect,
      fill = color_ecoli
    ) +
    
    geom_density(
      data = a.df[a.df[, "Species"] == "C.elegans", ],
      size = 0.6,
      alpha = 0.6,
      # key_glyph =  draw_key_rect,
      fill = color_celegans
    ) +
    

    
    
    labs(subtitle = a.subtitle, 
         caption = a.caption) +

 
    theme(
      panel.background  = element_blank(),
      plot.background = element_rect(fill = "white", color = "white"),
      
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        color = 'black',
        size = 12,
        angle = 0,
        margin = margin(0, 0, 0.0, 0)
      ),
      
      plot.subtitle = element_text(
        hjust = 0.5,
        face = "plain",
        color = 'black',
        size = 12,
        angle = 0,
        margin = margin(0, 0, 1, 0)
      ),
      
      plot.caption = element_text(
        hjust = -0.1,
        face = "italic",
        color = 'black',
        size = 7.5,
        angle = 0,
        margin = margin(5, 0,-4, 0)
      ),
      
      axis.text.x = element_text(
        face = "plain",
        color = 'black',
        size = 10,
        angle = 0,
        margin = margin(2, 0, 0, 0)
      ),
      
      axis.text.y = element_text(
        face = "plain",
        color = 'black',
        size = 10,
        angle = 0,
        margin = margin(0, 2, 0, 0)
      ),
      
      axis.title.y = element_text(
        face = "plain",
        color = 'black',
        size = 11,
        angle = 90,
        margin = margin(0, 4, 0, 0)
      ),
      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      axis.line = element_line(
        size = 0.3,
        linetype = "solid",
        color = "black"
      ),
      
      axis.ticks = element_line(
        size = 0.3,
        linetype = "solid",
        color = "black"
      ),
      
      axis.title.x = element_text(
        hjust = 0.5,
        face = "plain",
        color = 'black',
        size = 11,
        angle = 0,
        margin = margin(4, 0, 0, 0)
      ),
      
      
      legend.title = element_blank(),
      legend.background = element_rect(
        fill = "white",
        color = "white",
        size = 0.1),
      
      legend.key = element_rect(fill = "white", color = "white"),
      legend.key.height = unit(0, "cm"),
      legend.margin = margin(3, 0, 0, 0),
      legend.text.align = 0,
      legend.box.spacing = unit(0.0, "cm"),
      legend.box.margin = margin(0, 0, 0,-0.55, "cm"),
      legend.spacing.x = unit(0.0, "cm"),
      legend.position = "top",
      legend.text = element_text(
        face = "plain",
        color = "white",
        size = 9,
        angle = 0,
        margin = margin(0, 0, 0, 0))
 
      
    ) +
    
    guides(colour = guide_legend(override.aes = list(color = "white")))
      
      
      
    
}


# Generate plots.
plot_Prot_density <- f.p.density(Prot3,
                                 subtitle_Prot,
                                 caption_Prot)

plot_Prec_density <- f.p.density(Prec3,
                                 subtitle_Prec,
                                 caption_Prec)


# Export plots
f.export.plot(plot_Prot_density, "Prot_density", plot_width, plot_height, plot_res)
f.export.plot(plot_Prec_density, "Prec_density", plot_width, plot_height, plot_res)






# Plot simple boxplot ----------------------------

# This one is only combined with scatter and facet plots.
f.p.box_simple <- function(a.df) {
  ggplot(a.df,
         aes(x = Species,
             y = log2FC,
             fill = Species)) +
    
    scale_fill_manual(values = c(color_human, color_yeast, color_ecoli, color_celegans)) +
    
    labs(subtitle = "", caption = "") +
    
    coord_cartesian(ylim = c(FC_min, FC_max)) +
    scale_y_continuous(
      name = "",
      breaks = seq(from = round(digits = 0, FC_min), to = round(digits = 0, FC_max), by = 1 ),
      expand = expansion(mult = c(.02, .02))
    ) +
    
    theme(legend.position = "none") +
    theme_void() +
    theme(strip.text.x = element_blank()) +
    theme(text = element_text(color = "black")) +
    
    geom_boxplot(
      key_glyph = "blank",
      outlier.colour = "grey50",
      outlier.shape = 16,
      outlier.size = 1,
      outlier.alpha = 0.3,
      alpha = 0.7,
      notch = FALSE,
      fatten = 1,
      show.legend = FALSE
    )
}


# Generate but do not export plots.
plot_Prot_box_simple <- f.p.box_simple(Prot3)
plot_Prec_box_simple <- f.p.box_simple(Prec3)




# Plot_Scatter ------------------------------------------------------------

f.p.scatter <- function(a.df,
                        a.subtitle,
                        a.caption) {
  ggplot(a.df,
         aes(
           x = log2(a.df[, "ctr_mean"]),
           y = log2FC,
           color = Species,
           fill = Species
         )) +
    
    geom_hline(
      yintercept = median(a.df[a.df[, "Species"] == "Human", "log2FC"]),
      linetype = "dashed",
      color = color_human,
      alpha = 0.9
    ) +
    
    geom_hline(
      yintercept = median(a.df[a.df[, "Species"] == "Yeast", "log2FC"]),
      linetype = "dashed",
      color = color_yeast,
      alpha = 0.9
    ) +
    
    geom_hline(
      yintercept = median(a.df[a.df[, "Species"] == "E.coli", "log2FC"]),
      linetype = "dashed",
      color = color_ecoli,
      alpha = 0.9
    ) +
    
    geom_hline(
      yintercept = median(a.df[a.df[, "Species"] == "C.elegans", "log2FC"]),
      linetype = "dashed",
      color = color_celegans,
      alpha = 0.9
    ) +
    
    geom_point(shape = 16,
               size = 1.2,
               alpha = 0.20) +
    
    scale_color_manual(values = c(color_human, color_yeast, color_ecoli, color_celegans)) +
    scale_fill_manual(values = c(color_human, color_yeast, color_ecoli, color_celegans)) +
    
    labs(subtitle = a.subtitle,
         caption = a.caption) +
    
    coord_cartesian(ylim = c(FC_min, FC_max)) +
    scale_y_continuous(
      name = "log2(A:B)",
      breaks = seq(from = round(digits = 0, FC_min), to = round(digits = 0, FC_max), by = 1 ),
      expand = expansion(mult = c(.02, .02))
    ) +
    scale_x_continuous(name = "log2(B)",
                       breaks= pretty_breaks()) +
    
    theme(
      panel.background  = element_blank(),
      plot.background = element_rect(fill = "white", color = "white"),
      
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        color = 'black',
        size = 12,
        angle = 0,
        margin = margin(0, 0, 0.0, 0)
      ),
      
      plot.subtitle = element_text(
        hjust = 0.5,
        face = "plain",
        color = 'black',
        size = 12,
        angle = 0,
        margin = margin(0, 0, 1, 0)
      ),
      
      plot.caption = element_text(
        hjust = -0.1,
        face = "italic",
        color = 'black',
        size = 7.5,
        angle = 0,
        margin = margin(5, 0,-4, 0)
      ),
      
      axis.text.x = element_text(
        face = "plain",
        color = 'black',
        size = 10,
        angle = 0,
        margin = margin(2, 0, 0, 0)
      ),
      
      axis.text.y = element_text(
        face = "plain",
        color = 'black',
        size = 10,
        angle = 0,
        margin = margin(0, 2, 0, 0)
      ),
      
      axis.title.y = element_text(
        face = "plain",
        color = 'black',
        size = 11,
        angle = 90,
        margin = margin(0, 4, 0, 0)
      ),
      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      axis.line = element_line(
        size = 0.3,
        linetype = "solid",
        color = "black"
      ),
      
      axis.ticks = element_line(
        size = 0.3,
        linetype = "solid",
        color = "black"
      ),
      
      axis.title.x = element_text(
        hjust = 0.5,
        face = "plain",
        color = 'black',
        size = 11,
        angle = 0,
        margin = margin(4, 0, 0, 0)
      ),
      
      legend.title = element_blank(),
      legend.background = element_rect(
        fill = "transparent",
        color = NA,
        size = 0.1
      ),
      
      legend.key = element_rect(fill = "transparent", color = NA),
      legend.key.height = unit(0, "cm"),
      legend.margin = margin(3, 0, 0, 0),
      legend.text.align = 0,
      legend.box.spacing = unit(0.0, "cm"),
      legend.box.margin = margin(0, 0, 0,-0.55, "cm"),
      legend.spacing.x = unit(0.0, "cm"),
      legend.position = "top",
      legend.text = element_text(
        face = "plain",
        color = 'black',
        size = 9,
        angle = 0,
        margin = margin(0, 0, 0, 0)
      )
    ) +
    
    guides(colour = guide_legend(override.aes = list(alpha = 1)))
}


# Generate plots.
plot_Prot_scatter <-  f.p.scatter(Prot3,
                                  subtitle_Prot,
                                  caption_Prot)

plot_Prec_scatter <-  f.p.scatter(Prec3,
                                  subtitle_Prec,
                                  caption_Prec)

# Add boxplots on the side.
plot_Prot_scatter <-
  plot_grid(
    plot_Prot_scatter,
    plot_Prot_box_simple,
    labels = "",
    rel_widths = c(1, 0.25),
    align = "h",
    axis = "bt"
  )

plot_Prec_scatter <-
  plot_grid(
    plot_Prec_scatter,
    plot_Prec_box_simple,
    labels = "",
    rel_widths = c(1, 0.25),
    align = "h",
    axis = "bt"
  )


# Export plots.
f.export.plot(plot_Prot_scatter, "Prot_scatter", plot_width, plot_height, plot_res)
f.export.plot(plot_Prec_scatter, "Prec_scatter", plot_width, plot_height, plot_res)





# Plot facet --------------------------------------------------------------

f.p.facet <- function(a.df,
                      a.subtitle,
                      a.caption) {
  ggplot(a.df,
         aes(
           x = log2(a.df[, "ctr_mean"]),
           y = log2FC,
           color = Species,
           fill = Species
         )) +
    
    geom_hline(
      yintercept = median(a.df[a.df[, "Species"] == "E.coli", "log2FC"]),
      linetype = "dashed",
      color = color_ecoli,
      alpha = 0.7
    ) +

    geom_hline(
      yintercept = median(a.df[a.df[, "Species"] == "Human", "log2FC"]),
      linetype = "dashed",
      color = color_human,
      alpha = 0.7
    ) +

    geom_hline(
      yintercept = median(a.df[a.df[, "Species"] == "Yeast", "log2FC"]),
      linetype = "dashed",
      color = color_yeast,
      alpha = 0.7
    ) +
    
    geom_hline(
      yintercept = median(a.df[a.df[, "Species"] == "C.elegans", "log2FC"]),
      linetype = "dashed",
      color = color_celegans,
      alpha = 0.7
    ) +
    
    
    geom_point(shape = 16,
               size = 1.2,
               alpha = 0.3) +
    scale_color_manual(values = c(color_human, color_yeast, color_ecoli, color_celegans)) +
    scale_fill_manual(values = c(color_human, color_yeast, color_ecoli, color_celegans)) +
    
    labs(subtitle = a.subtitle,
         caption = a.caption) +
    
    theme(
      panel.background  = element_blank(),
      plot.background = element_rect(fill = "white", color = "white"),
      
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        color = 'black',
        size = 12,
        angle = 0,
        margin = margin(0, 0, 0.0, 0)
      ),
      
      plot.subtitle = element_text(
        hjust = 0.5,
        face = "plain",
        color = 'black',
        size = 12,
        angle = 0,
        margin = margin(0, 0, 1, 0)
      ),
      
      plot.caption = element_text(
        hjust = -0.1,
        face = "italic",
        color = 'black',
        size = 7.5,
        angle = 0,
        margin = margin(5, 0,-4, 0)
      ),
      
      axis.text.x = element_text(
        face = "plain",
        color = 'black',
        size = 10,
        angle = 0,
        margin = margin(2, 0, 0, 0)
      ),
      
      axis.text.y = element_text(
        face = "plain",
        color = 'black',
        size = 10,
        angle = 0,
        margin = margin(0, 2, 0, 0)
      ),
      
      axis.title.y = element_text(
        face = "plain",
        color = 'black',
        size = 11,
        angle = 90,
        margin = margin(0, 4, 0, 0)
      ),
      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      axis.line = element_line(
        size = 0.3,
        linetype = "solid",
        color = "black"
      ),
      
      axis.ticks = element_line(
        size = 0.3,
        linetype = "solid",
        color = "black"
      ),
      
      axis.title.x = element_text(
        hjust = 0.5,
        face = "plain",
        color = 'black',
        size = 11,
        angle = 0,
        margin = margin(4, 0, 0, 0)
      ),
      
      legend.title = element_blank(),
      legend.background = element_rect(
        fill = "transparent",
        color = NA,
        size = 0.1
      ),
      
      legend.key = element_rect(fill = "transparent", color = NA),
      legend.key.height = unit(0, "cm"),
      legend.margin = margin(3, 0, 0, 0),
      legend.text.align = 0,
      legend.box.spacing = unit(0.0, "cm"),
      legend.box.margin = margin(0, 0, 0,-0.55, "cm"),
      legend.spacing.x = unit(0.0, "cm"),
      legend.position = "top",
      legend.text = element_text(
        face = "plain",
        color = 'black',
        size = 9,
        angle = 0,
        margin = margin(0, 0, 0, 0)
      )
    ) +
    coord_cartesian(ylim = c(FC_min, FC_max)) +
    scale_y_continuous(
      name = "log2(A:B)",
      breaks = seq(from = round(digits = 0, FC_min), to = round(digits = 0, FC_max), by = 1 ),
      expand = expansion(mult = c(.02, .02))
    ) +
    scale_x_continuous(name = "log2(B)",
                       expand = expansion(mult = c(.02, .02)),
                       breaks = pretty_breaks(),
                       )  +
    facet_grid(. ~ Species) +
    theme(strip.text.x = element_blank()) +
    guides(colour = guide_legend(override.aes = list(alpha = 1)))
}


# Generate Plots
plot_Prot_facet_ori <-  f.p.facet(Prot3,
                              subtitle_Prot,
                              caption_Prot)

plot_Prec_facet_ori <-  f.p.facet(Prec3,
                              subtitle_Prec,
                              caption_Prec)


# Add boxplot on the side
plot_Prot_facet <-
  plot_grid(
    plot_Prot_facet_ori,
    plot_Prot_box_simple,
    labels = "",
    rel_widths = c(1, 0.25),
    align = "h",
    axis = "bt"
  )

plot_Prec_facet <-
  plot_grid(
    plot_Prec_facet_ori,
    plot_Prec_box_simple,
    labels = "",
    rel_widths = c(1, 0.25),
    align = "h",
    axis = "bt"
  )


# Export combined plots
f.export.plot(plot_Prot_facet, "Prot_facet", plot_width, plot_height, plot_res)
f.export.plot(plot_Prec_facet, "Prec_facet", plot_width, plot_height, plot_res)





# Plot Volcano ------------------------------------------------------------

f.p.volcano <- function(a.df,
                        a.subtitle,
                        a.caption)   {
  ggplot(a.df,
         aes(
           x = log2FC,
           y = -log10(p_adj),
           color = Species,
           fill = Species
         )) +
    

    geom_hline(
      yintercept = -log10(alpha_limma),
      linetype = "dotted",
      color = "grey60",
      alpha = 0.8
    ) +
    
    geom_vline(
      xintercept = c(-limit_FC,+limit_FC),
      linetype = c("dotted", "dotted"),
      color = c("grey60", "grey60"),
      alpha = c(0.8, 0.8)
    ) +
    
    geom_point(shape = 16,
               size = 1.2,
               alpha = 0.3) +
    
    scale_color_manual(values = c(color_human, color_yeast, color_ecoli, color_celegans)) +
    scale_fill_manual(values = c(color_human, color_yeast, color_ecoli, color_celegans)) +
    
    labs(subtitle = a.subtitle,
         caption = a.caption) +
    
    theme(
      panel.background  = element_blank(),
      plot.background = element_rect(fill = "white", color = "white"),
      
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        color = 'black',
        size = 12,
        angle = 0,
        margin = margin(0, 0, 0.0, 0)
      ),
      
      plot.subtitle = element_text(
        hjust = 0.5,
        face = "plain",
        color = 'black',
        size = 12,
        angle = 0,
        margin = margin(0, 0, 1, 0)
      ),
      
      plot.caption = element_text(
        hjust = -0.1,
        face = "italic",
        color = 'black',
        size = 7.5,
        angle = 0,
        margin = margin(5, 0,-4, 0)
      ),
      
      axis.text.x = element_text(
        face = "plain",
        color = 'black',
        size = 10,
        angle = 0,
        margin = margin(2, 0, 0, 0)
      ),
      
      axis.text.y = element_text(
        face = "plain",
        color = 'black',
        size = 10,
        angle = 0,
        margin = margin(0, 2, 0, 0)
      ),
      
      axis.title.y = element_text(
        face = "plain",
        color = 'black',
        size = 11,
        angle = 90,
        margin = margin(0, 4, 0, 0)
      ),
      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      axis.line = element_line(
        size = 0.3,
        linetype = "solid",
        color = "black"
      ),
      
      axis.ticks = element_line(
        size = 0.3,
        linetype = "solid",
        color = "black"
      ),
      
      axis.title.x = element_text(
        hjust = 0.5,
        face = "plain",
        color = 'black',
        size = 11,
        angle = 0,
        margin = margin(4, 0, 0, 0)
      ),
      
      legend.title = element_blank(),
      legend.background = element_rect(
        fill = "transparent",
        color = NA,
        size = 0.1
      ),
      
      legend.key = element_rect(fill = "transparent", color = NA),
      legend.key.height = unit(0, "cm"),
      legend.margin = margin(3, 0, 5, 0),
      legend.text.align = 0,
      legend.box.spacing = unit(0.0, "cm"),
      legend.box.margin = margin(0, 0, 0,-0.55, "cm"),
      legend.spacing.x = unit(0.0, "cm"),
      legend.position = "top",
      legend.text = element_text(
        face = "plain",
        color = 'black',
        size = 9,
        angle = 0,
        margin = margin(0, 0, 0, 0)
      )
    ) +
    
    coord_cartesian(xlim = c(FC_min, FC_max),
                    ylim = c(0,10)) +
    scale_x_continuous(
      name = "log2(A:B)",
      breaks = seq(from = round(digits = 0, FC_min), to = round(digits = 0, FC_max), by = 1 ),
      expand = expansion(mult = c(.02, .02))
    ) +
    scale_y_continuous(name = "-log10(p_adj)") +
    guides(colour = guide_legend(override.aes = list(alpha = 1)))
}


# Generate plot.
plot_Prot_volcano <-  f.p.volcano(Prot3,
                                  subtitle_Prot,
                                  caption_Prot)

# Export plot.
f.export.plot(plot_Prot_volcano, "Prot_volcano", plot_width, plot_height, plot_res)





##### The End