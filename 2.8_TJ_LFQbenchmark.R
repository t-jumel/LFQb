# Overview ------------------------------------------------------------

### Literature:
# Demichev, Vadim, et al. "DIA-NN: neural networks and interference correction enable deep proteome coverage in high throughput." Nature methods 17.1 (2020): 41-44.
# Kuharev, Jörg, et al. "In-depth evaluation of software tools for data-independent acquisition based label-free quantification." Proteomics 15.18 (2015): 3140-3151.
# Ritchie, Matthew E., et al. "limma powers differential expression analyses for RNA-sequencing and microarray studies." Nucleic acids research 43.7 (2015): e47-e47.


### Disclaimer:
# Hi, my name is tobias and I hope you can make use
# of this script. Please excuse if the code is not elegant
# as I am a beginner with R. I am thankful for suggestions.
# This script was inspired by the the LFQbench paper
# (Kuharev et al.), and serves as alternative for the related package.


### Contact:
# Tobias Jumel
# jumel@mpi-cbg.de
# Shevchenko Lab, MPI-CBG
# 14.06.2022
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

# Select a folder location (variable "folder_input").
# This folder must contain exactly one protein group
# and precursor .tsv matrix in a DIA-NN style layout (column names)
# with "pg_matrix" and "pr_matrix" in the respective filenames.
# "first-run" matrices need to be be moved to subfolders 
# or the import function edited.

# Check package installations, especially bioconductor ones.

# Specify all parameters under "variables" section, in most cases
# this applies just to filter settings and column indices.
# Then execute and find the output folder (csv files, plots)
# within the specified input folder.




# Introduction -------------------------------------------------------------

### General
# This script analyses benchmarks of peptide and protein group
# quantification by bottom-up proteomics with LC-MS and DIA-NN.


# The benchmarks use 2 samples, each consisting of multiple species
# mixed together in different percentages. 
# Parts of this script are hard-coded matching the example below!
# classical LFQbenchmark sample A: 65% Human, 30% Yeast, 5% E.coli
# classical LFQbenchmark sample B: 65% Human, 15% Yeast, 20% E.coli


# When comparing peptide and protein quantities from sample B to A,
# one expects human proteins to be stable, yeast 2x up, and E.coli
# 4x down regulated. The benchmarks reveal to which extent 
# measured log2 fold-changes deviate from expected log2FC,
# including e.g. global over or underestimation and identification errors.


# Note: This script should directly work on 2-species benchmarks,
# but the interpretation of up- and downregulation is hard-coded,
# (function f.LFQbenchmark.limma.interpretation)
# Try use sample setups close to the introduction example,
# otherwise the deFDR and other summary stats might produce errors.



### Limma and LFQ benchmarks
# This script additionally uses differential expression analysis,
# meaning a statistical test for each protein whether it is upregulated,
# downregulated, or stable. Matching these test results against expectations
# derived form the defined sample mixtures in a "confusion matrix" style 
# provides the most useful summary statistics (TP and deFDR, see below).
# These among others can be used to evaluate the performance and reliability
# of multiple result sets in comparison (highly useful to optimize DIA-NN settings).

# (e.g. A yeast protein group measured as upregulated is a true-positive (TP),
# a human protein group measured as downregulated a false-positive, etc...)


### Interpretation
# A false-discovery rate (deFDR) is calculated from those positives
# and negatives (confusion matrix layout, deFDR = TP/(TP+FP)).
# It can indicate problems with the reliability of a workflow
# if exceeding the threshold of the statistical limma analysis. 
# (alpha_limma variable, typically 0.01 or 1 %)


### Interpretation Tips
# E.g. if alpha_limma = 0.01 and the resulting deFDR is 1.9% this indicates
# considerable issues within the data,
# such as false identifications, erroneous normalization,
# or global quantification errors such as fold-change overestimation.

# If multiple datasets result in similar, acceptable deFDR values (often <1%),
# and are otherwise similar, the one with the highest TP count is "the best".
# The TP stat is more useful than the number of identified and quantified proteins,
# as it includes not only sensitivity aspects but also quantitative precision
# and to some degree quantitative accuracy issues like ratio compression 
# (fold-change underestimation, somewhat common in TOF data).

# The plots (especially facet and density) can visualize issues related to
# excess false-positives and dataset-wide tendencies for over/underestimation.
# Check for individual points far away from the expected fold-changes 
# especially if those errors seem to cluster together at log2 fold-changes
# matching a different species. (Potentially false signal-to-peptide matches, 
# inflated precursor or protein group FDR, etc.)

# Also check for proper normalization and over- or underestimation of
# fold-changes on density plots, especially in 2-species benchmarks.




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



# Pre-typed to install bioconductor packages.


# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("limma")
# BiocManager::install("Biobase")

library(limma)
library(Biobase)




# Variables ---------------------------------------------------------------

# Explanation below at section end.


# Set working directory (just pre-typed, might not be used).
# setwd()

# Script version and name used in output filenames.
script_name <- "TJ_LFQbenchmark"
script_version <- "2.8"

# Recommend use of "path copy copy".
# See "To run this Script" section for input requirements.
folder_input <- "//fileserver/project_jumel/All_LFQ_benchmarks/Reanalysis_external/ReAnalysis_LFQbench_PXD028735/diaPASEF_data/diaPASEF03"



# Naming the 2 different benchmark sample mixtures.
# This is just for record-keeping and preventing errors 
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


# p_adj cut-off for diff. expr. analysis by limma.
alpha_limma <- 0.01


# Expected fold-changes by sample mixtures.
expFC_human <- 0
expFC_yeast <- +1
expFC_ecoli <- -2


# Variables for plots.

# Colors
color_human <- "#199d76"
color_yeast <- "#d85f02"
color_ecoli <- "#7570b2"

# These plot limits are set a bit higher (0.2) than 
# needed to prevent value clipping.
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
# in both cases for experimental and control condition.

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
# considered as "NOT sign. up- or downregulated",
# no matter the limma result (p_adj).
# This value might be set ato be between 0.5-1.0 .

# alpha_limma is the stat. significance cutoff 
# for the p_adj of limma. Always compare the resulting deFDR to this value,
# an excess false-positive rate indicates potential errors in the data.

# expFC_ecoli and similar refer the species-specific 
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
    folder_input, "/",script_version,"_",script_name,
    "_", rep_min,"of", rep_max,
    "_CV", limit_CV,
    "_FC",limit_FC,
    "_a",alpha_limma)

dir.create(folder_output)



# Collect variables in dataframe for export.
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
    color_human,
    color_yeast,
    color_ecoli,
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
  "color_human",
  "color_yeast",
  "color_ecoli",
  "folder_output"
)

# Export variables as separate .csv file.
write.csv(variables,
          row.names = TRUE,
          paste0(folder_output, "/", "variables_",
                 basename(dirname(folder_output)), ".csv"))



# Functions ---------------------------------------------------------------

# Recognize user defined functions (f.xyz.abc) and arguments (a.xyz.abc)
# Plot functions are at the script end (f.p.scatter) together with
# the respective plot export functions.


# Import DIA-NN matrix from specified folder into dataframe.
# Uses partial filename match, e.g."pg_matrix", "pr_matrix".
# requires e.g. first-pass matrices to be moved to subfolders.
# e.g. Prot0 <- f.import.matrix(folder_input, "pg_matrix", "\t")
f.import.matrix <- function(a.folder,
                            a.matrix,
                            a.sep) {
  read.csv(sep = a.sep,
           list.files(a.folder,
                      pattern = a.matrix,
                      full.names = TRUE))
}



# In DIA-NN matrices, rename columns containing quantitative values
# from filepath to short (stripped) sample name.
# Remove file suffix and filepath.
# e.g. Prot1 <- f.strip.sample.column.name(Prot1, col_ctr_Prot)
f.strip.sample.column.name <- function(a.df, a.columns) {
  names(a.df)[a.columns] <-
    sub('.*\\.', '', gsub(".mzML|.raw|.dia", "", colnames(a.df[a.columns])))
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
# This contains also annotating entries with respective species 
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
  
  a.df[, "Species"][grepl('_ECOLI', a.df[, "Protein.Names"]) &
                      !grepl('_HUMAN', a.df[, "Protein.Names"]) &
                      !grepl('_YEAST', a.df[, "Protein.Names"])] <-
    'E.coli'
  
  a.df[, "Species"][grepl('_HUMAN', a.df[, "Protein.Names"]) &
                      !grepl('_YEAST', a.df[, "Protein.Names"]) &
                      !grepl('_ECOLI', a.df[, "Protein.Names"])] <-
    'Human'
  
  a.df[, "Species"][grepl('_YEAST', a.df[, "Protein.Names"]) &
                      !grepl('_HUMAN', a.df[, "Protein.Names"]) &
                      !grepl('_ECOLI', a.df[, "Protein.Names"])] <-
    'Yeast'
  
  a.df <- a.df[complete.cases(a.df[, "Species"]), ]
  
  
  
  # Add column listing expected log2FC row-wise
  # to simplify calculations of summary_stats around the log2FC.
  a.df[, "expected_log2FC"] <-
    ifelse(a.df[, "Species"] == "Human",
           expFC_human,
           ifelse(a.df[, "Species"] == "Yeast",
                  expFC_yeast,
                  expFC_ecoli))
  
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
  grps = c(rep('exp', length(a.col_exp)), rep('con', length(a.col_ctr))) # conditions of samples in assay
  reps = c(1:length(a.col_exp), 1:length(a.col_ctr))  # replicates
  
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
  limma.rep <- factor(reps)
  contrast.matrix <- model.matrix(~  0 + limma.cond + limma.rep)
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
  
  
  # Add p_adj (limma result) to main protein group dataframe.
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




# Interpret the P_adj from limma differential expression analysis
# by comparing results against expectations from the sample mixtures.
# (yeast protein groups found as upregulated is annotated as true-positive,
# human protein groups found as upregulated is annotated as false-positive etc.)

# !!! hard-coded to the example from the introduction.
# While this tolerates different extents of fold-changes,
# the direction of up and downregulation are fixed to the mixtures
# as describe din the introduction.

# If you deviate from (B to A: human-stable, yeast up, E.coli down);
# this function might produce nonsense errors. 2-species mixtures with the same
# "direction" of regulation should work without adjustment.

# The resulting true and false - positives and negatives can be seen as
# "confusion matrix" summary stats and are of highest value to evaluate results.
# e.g. Prot3 <- f.LFQbenchmark.limma.interpretation(Prot3, alpha_limma, limit_FC)
f.LFQbenchmark.limma.interpretation <-
  function(a.df, a.alpha_limma, a.limit_FC) {
    
    a.df[(a.df[, "Species"] == "Human"
          & (a.df[, "p_adj"] >= a.alpha_limma
             | abs(a.df[, "log2FC"]) <= a.limit_FC)), "DE_result"] <-
      "true-negative"
    
    a.df[(a.df[, "Species"] == "Human"
          & a.df[, "p_adj"] < a.alpha_limma
          & abs(a.df[, "log2FC"]) > a.limit_FC), "DE_result"] <-
      "false-positive"
    
    
    a.df[(a.df[, "Species"] == "Yeast"
          & a.df[, "log2FC"] > +a.limit_FC
          & a.df[, "p_adj"] < a.alpha_limma), "DE_result"] <-
      "true-positive"
    
    a.df[(a.df[, "Species"] == "Yeast"
          & a.df[, "log2FC"] < -a.limit_FC
          & a.df[, "p_adj"] < a.alpha_limma), "DE_result"] <-
      "false-positive"
    
    a.df[(a.df[, "Species"] == "Yeast"
          & (a.df[, "p_adj"] >= a.alpha_limma
             | abs(a.df[, "log2FC"]) <= a.limit_FC)), "DE_result"] <-
      "false-negative"
    
    
    a.df[(a.df[, "Species"] == "E.coli"
          & a.df[, "log2FC"] < -a.limit_FC
          & a.df[, "p_adj"] < a.alpha_limma), "DE_result"] <-
      "true-positive"
    
    a.df[(a.df[, "Species"] == "E.coli"
          & a.df[, "log2FC"] > +a.limit_FC
          & a.df[, "p_adj"] < a.alpha_limma), "DE_result"] <-
      "false-positive"
    
    a.df[(a.df[, "Species"] == "E.coli"
          & (a.df[, "p_adj"] >= a.alpha_limma
             | abs(a.df[, "log2FC"]) <= a.limit_FC)), "DE_result"] <-
      "false-negative"
    
    return(a.df)
    
  }




# For LFQbenchmark data, calculate skewness and kurtosis for both
# precursor and protein group level, collect in separate dataframe.
# While not the most useful stats, high values can hint at potential problems
# of the quantification like biases. Use with caution and just as an indicator.
# e.g. skew_kurt <-  f.LFQbenchmark.skewness(Prec3, Prot3, "log2FC")
f.LFQbenchmark.skewness <- function(a.df01, a.df02, a.column) {
  
  skew_kurt <- data.frame(matrix(ncol = 0, nrow = 6))
  skew_kurt[, "Species"] <- c("E.coli", "Human", "Yeast")
  skew_kurt[, "Group"] <- c(rep("precursor", 3), 
                            rep("protein group", 3))
  
  skew_kurt[, "Skewness"] <- c(
    round(digits = 2, skewness(a.df01[which(a.df01[, "Species"] == "E.coli"), a.column])),
    round(digits = 2, skewness(a.df01[which(a.df01[, "Species"] == "Human"), a.column])),
    round(digits = 2, skewness(a.df01[which(a.df01[, "Species"] == "Yeast"), a.column])),
    round(digits = 2, skewness(a.df02[which(a.df02[, "Species"] == "E.coli"), a.column])),
    round(digits = 2, skewness(a.df02[which(a.df02[, "Species"] == "Human"), a.column])),
    round(digits = 2, skewness(a.df02[which(a.df02[, "Species"] == "Yeast"), a.column]))
  )
  
  skew_kurt[, "Kurtosis"] <- c(
    round(digits = 2, kurtosis(a.df01[which(a.df01[, "Species"] == "E.coli"), a.column])),
    round(digits = 2, kurtosis(a.df01[which(a.df01[, "Species"] == "Human"), a.column])),
    round(digits = 2, kurtosis(a.df01[which(a.df01[, "Species"] == "Yeast"), a.column])),
    round(digits = 2, kurtosis(a.df02[which(a.df02[, "Species"] == "E.coli"), a.column])),
    round(digits = 2, kurtosis(a.df02[which(a.df02[, "Species"] == "Human"), a.column])),
    round(digits = 2, kurtosis(a.df02[which(a.df02[, "Species"] == "Yeast"), a.column]))
  )
  
  return(skew_kurt)
}




# For processed LFQbenchmark data, collect summary statistics
# into new dataframe for export and plot subtitles.
# They describe sensitivity, precision, accuracy, reliability etc..
# Next to plots they are the core result of this script,
# mainly use true-positives (TP) and deFDR to evaluate results.
f.LFQbenchmark.summary.stats <-
  function(a.df01, a.df02, a.df03, a.df04) {
    summary_stats <- data.frame(matrix(ncol = 0, nrow = 1))
    
    
    # Number of confidently detected and quantified protein groups,
    # after filtering for missingness and CV, respectively.
    summary_stats[, "Prot_ID"] <- nrow(a.df01)
    summary_stats[, "Prot_Quant"] <- nrow(a.df02)

    
    # deFDR (differential expression FDR) ("confusion matrix)
    # Useful to spot different prevalence of aberrant quantifications.
    # Does not recognize mismatch errors within one species.
    # Especially useful to spot errors which survive typical filtering
    # like missingness or CV such as false signal to peptide assignments.
    # If alpha_limma is set to 0.01, the deFDR should be around 1%
    # 2% would already indicate severe problems in the data.
    summary_stats[, "deFDR"] <-
      round(digits = 2,
            100 * length(which(a.df02[, "DE_result"] == "false-positive")) /
              (length(which(
                a.df02[, "DE_result"] == "false-positive"
              ))
              + length(which(
                a.df02[, "DE_result"] == "true-positive"
              ))))
    
    
    # Confusion matrix of differential expression results.
    # Secondary stats are continued below.
    summary_stats[, "TP"] <-
      length(which(a.df02[, "DE_result"] == "true-positive"))
    
    
    # Accuracy is average distance of measured to expected log2FC.
    summary_stats[, "Prot_Accuracy"] <-
      round(digits = 2,
            mean(abs(a.df02[, "log2FC"] - a.df02[, "expected_log2FC"]))
            )
    
    # Trueness is distance between measurement medians 
    # and respective expected fold-changes.
    # (cumulative as shifts are important to spot)
    # Several conditions can lead to a high value
    # e.g. ratio compression or erroneous normalization or even both together.
    summary_stats[, "Prot_Trueness"] <-
      round(
        digits = 2,
        sum(na.rm = TRUE,
        c(abs(median(a.df02[a.df02[, "Species"] == "E.coli" , "log2FC"]) - expFC_ecoli),
          abs(median(a.df02[a.df02[, "Species"] == "Human"  , "log2FC"]) - expFC_human),
          abs(median(a.df02[a.df02[, "Species"] == "Yeast"  , "log2FC"]) - expFC_yeast))))
    
    
    
    # Dispersion is the average distance of individual measurements around the respective median.
    # Can be used to evaluate quantification in case the difference between median measurements
    # and expectations is high, which means differential expression analysis is not meaningful.
    summary_stats[, "Prot_Dispersion"] <-
    round(digits = 2,
          mean(na.rm = TRUE,
               c(
                 abs(a.df02[a.df02[, "Species"] == "Human"   , "log2FC"] - median(a.df02[a.df02[, "Species"] == "Human"  , "log2FC"])),
                 abs(a.df02[a.df02[, "Species"] == "Yeast"   , "log2FC"] - median(a.df02[a.df02[, "Species"] == "Yeast"  , "log2FC"])),
                 abs(a.df02[a.df02[, "Species"] == "E.coli"  , "log2FC"] - median(a.df02[a.df02[, "Species"] == "E.coli" , "log2FC"]))
               )))

    
    # Coefficient of variation to describe precision.
    # For ca. 3 replicates aim at average and median at or below 5%
    # after filtering for <20%.
    summary_stats[, "Prot_CV_Mean"]   <-
      round(digits = 2,  mean(c(a.df02[, "ctr_CV"], a.df02[, "exp_CV"])))
    summary_stats[, "Prot_CV_Median"] <-
      round(digits = 2, median(c(a.df02[, "ctr_CV"], a.df02[, "exp_CV"])))
    
    

    
    # Other than TP, continued confusion matrix of 
    # differential expression result interpretation.
    summary_stats[, "FP"] <-
      length(which(a.df02[, "DE_result"] == "false-positive"))
    summary_stats[, "TN"] <-
      length(which(a.df02[, "DE_result"] == "true-negative"))
    summary_stats[, "FN"] <-
      length(which(a.df02[, "DE_result"] == "false-negative"))
    
    summary_stats[, "Sensitivity"] <-
      round(digits = 2, 100 * length(which(a.df02[, "DE_result"] == "true-positive")) /
              (length(which(
                a.df02[, "DE_result"] == "true-positive"
              )) + length(which(
                a.df02[, "DE_result"] == "false-negative"
              ))))
    
    summary_stats[, "Specificity"] <-
      round(digits = 2, 100 * length(which(a.df02[, "DE_result"] == "true-negative")) /
              (length(which(
                a.df02[, "DE_result"] == "true-negative"
              )) + length(which(
                a.df02[, "DE_result"] == "false-positive"
              ))))
    
    
    # Skewness, just the most important values of that section collected 
    # here, all values are in a separate skew_kurt dataframe 
    # and exported table.
    summary_stats[, "Skewness_Prot_E.coli"] <-
      round(digits = 2, skewness(a.df02[which(a.df02[, "Species"] == "E.coli"), "log2FC"]))
    summary_stats[, "Skewness_Prot_Yeast"]  <-
      round(digits = 2, skewness(a.df02[which(a.df02[, "Species"] == "Yeast"), "log2FC"]))
    
    
    
    
    # Stats partially repeated on precursor level as before.
    summary_stats[, "Prec_ID"] <- nrow(a.df03)
    summary_stats[, "Prec_Quant"] <- nrow(a.df04)

    

    # Accuracy is average distance of measured to expected log2FC.
    summary_stats[, "Prec_Accuracy"] <-
      round(digits = 2,
            mean(abs(a.df04[, "log2FC"] - a.df04[, "expected_log2FC"]))
      )
    
    # Trueness is distance between measurement medians 
    # and respective expected fold-changes (cumulative as shifts are important to spot)
    # Several conditions can lead to a high value
    # e.g. ratio compression or erroneous normalization or even both together.
    summary_stats[, "Prec_Trueness"] <-
      round(
        digits = 2,
        sum(na.rm = TRUE,
            c(abs(median(a.df04[a.df04[, "Species"] == "E.coli"  , "log2FC"]) - expFC_ecoli),
              abs(median(a.df04[a.df04[, "Species"] == "Human"  , "log2FC"]) - expFC_human),
              abs(median(a.df04[a.df04[, "Species"] == "Yeast"  , "log2FC"]) - expFC_yeast))))
    
    
    
    # Dispersion is the average distance of individual measurements around the respective median.
    # Can be used to evaluate quantification in case the difference between median measurements
    # and expectations is high, which means differential expression analysis is not meaningful.
    summary_stats[, "Prec_Dispersion"] <-
    round(digits = 2,
          mean(na.rm = TRUE,
               c(
                 abs(a.df04[a.df04[, "Species"] == "Human"  , "log2FC"] - median(a.df04[a.df04[, "Species"] == "Human"  , "log2FC"])),
                 abs(a.df04[a.df04[, "Species"] == "Yeast"  , "log2FC"] - median(a.df04[a.df04[, "Species"] == "Yeast"  , "log2FC"])),
                 abs(a.df04[a.df04[, "Species"] == "E.coli"  , "log2FC"] - median(a.df04[a.df04[, "Species"] == "E.coli"  , "log2FC"]))
               )))
    
    
    summary_stats[, "Prec_CV_Mean"]   <-  round(digits = 2,
                                                mean(c(a.df04[, "ctr_CV"], a.df04[, "exp_CV"])))
    summary_stats[, "Prec_CV_Median"] <-  round(digits = 2,
                                                median(c(a.df04[, "ctr_CV"], a.df04[, "exp_CV"])))
    
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
Prot2 <- subset(Prot1, ctr_count >= rep_min & exp_count >= rep_min)


# Retain entries matching defined precision cutoff.
Prot3 <- subset(Prot2, ctr_CV <  limit_CV & exp_CV <  limit_CV)


# Differential expression analysis and match of results
# against expectations from sample mixtures.
Prot3 %<>%
  f.limma.p_adj(col_ctr_Prot, col_exp_Prot) %>%
  f.LFQbenchmark.limma.interpretation(alpha_limma, limit_FC)



# Repeat data processing except differential expression
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
Prec2 <- subset(Prec1, ctr_count >= rep_min & exp_count >= rep_min)


# Retain entries matching defined precision cutoff.
Prec3 <- subset(Prec2, ctr_CV <  limit_CV & exp_CV <  limit_CV)







# Summary Stats -----------------------------------


# Create extra dataframe containing skewness ans kurtosis values.
# A subselection is also present in the main summary_stat dataframe and export.
skew_kurt <-  f.LFQbenchmark.skewness(Prec3, Prot3, "log2FC")

# Export skew_kurt.
write.csv(
  skew_kurt,
  paste0(
    folder_output,
    "/",
    "skewness_kurtosis",
    basename(dirname(folder_output)),
    ".csv"
  ),
  row.names = FALSE
)



# Create dataframe of LFQbenchmark summary stats (main result)
summary_stats <- f.LFQbenchmark.summary.stats(Prot2, Prot3, Prec2, Prec3)

# Export summary stats.
write.csv(
  summary_stats,
  paste0(
    folder_output,
    "/",
    "summary_stats_",
    basename(dirname(folder_output)),
    ".csv"
  ),
  row.names = FALSE
)










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
  " true-positives, ",
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
  factor(Prot3$Species, levels = c("Human", "Yeast", "E.coli"))
Prec3$Species <-
  factor(Prec3$Species, levels = c("Human", "Yeast", "E.coli"))





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
    xlim = c(1.05*min(pca_df[,"x"]),1.05*max(pca_df[,"x"])),
    ylim = c(1.05*min(pca_df[,"y"]),1.05*max(pca_df[,"y"])))+
  
  geom_point(size = 2.5, alpha = 0.5) +
  
  geom_text(size = 2.5, 
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



# Plot Skewness, Kurtosis -----------------------------------------------

f.p.skewness <- function(a.df,
                         a.subtitle,
                         a.caption) {
  ggplot(a.df,
         aes(
           x = Skewness,
           y = Kurtosis,
           color = Species,
           fill = Species,
           shape = Group
         )) +
    
    geom_point(size = 2.0, alpha = 0.7) +
    
    scale_color_manual(values = c(color_human, color_yeast, color_ecoli)) +
    scale_fill_manual(values = c(color_human, color_yeast, color_ecoli)) +
    
    labs(subtitle = a.subtitle,
         caption = a.caption) +
    
    coord_cartesian(xlim = c(-max(abs(skew_kurt[, "Skewness"])),+max(abs(skew_kurt[, "Skewness"])))) +
    
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
plot_skewness <-  f.p.skewness(skew_kurt,
                               subtitle_Prot,
                               caption_Prot)

# Export plots.
f.export.plot(plot_skewness, "Prot_skewness", plot_width, plot_height, plot_res)











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
    
    scale_fill_manual(values = c(color_human, color_yeast, color_ecoli)) +
    
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


    
    geom_density(
      data = a.df[a.df[, "Species"] == "Human", ],
      size = 0.6,
      alpha = 0.6,
      key_glyph =  draw_key_point,
      fill = color_human
    ) +
    
    geom_density(
      data = a.df[a.df[, "Species"] == "Yeast", ],
      size = 0.6,
      alpha = 0.6,
      key_glyph =  draw_key_point,
      fill = color_yeast
      ) +
    
    geom_density(
      data = a.df[a.df[, "Species"] == "E.coli", ],
      size = 0.6,
      alpha = 0.6,
      key_glyph =  draw_key_point,
      fill = color_ecoli
    ) +
    
    labs(subtitle = a.subtitle, caption = a.caption) +

 
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
      legend.margin = margin(2, 0, 5, 0),
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
    )
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
    
    scale_fill_manual(values = c(color_human, color_yeast, color_ecoli)) +
    
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
    
    geom_point(shape = 16,
               size = 1.2,
               alpha = 0.20) +
    
    scale_color_manual(values = c(color_human, color_yeast, color_ecoli)) +
    scale_fill_manual(values = c(color_human, color_yeast, color_ecoli)) +
    
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
    
    
    geom_point(shape = 16,
               size = 1.2,
               alpha = 0.3) +
    scale_color_manual(values = c(color_human, color_yeast, color_ecoli)) +
    scale_fill_manual(values = c(color_human, color_yeast, color_ecoli)) +
    
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
    
    scale_color_manual(values = c(color_human, color_yeast, color_ecoli)) +
    scale_fill_manual(values = c(color_human, color_yeast, color_ecoli)) +
    
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