# Overview ------------------------------------------------------------



### General
# This script analyses benchmarks of peptide precursor and protein group
# (mostly label-free) quantification by bottom-up proteomics 
# with DIA LC-MS and DIA-NN.



### Literature:
# Demichev, Vadim, et al. "DIA-NN: neural networks and interference correction enable deep proteome coverage in high throughput." Nature methods 17.1 (2020): 41-44.
# Kuharev, Jörg, et al. "In-depth evaluation of software tools for data-independent acquisition based label-free quantification." Proteomics 15.18 (2015): 3140-3151.
# Navarro, Pedro, et al. "A multicenter study benchmarks software tools for label-free proteome quantification." Nature biotechnology 34.11 (2016): 1130-1136.
# Ritchie, Matthew E., et al. "limma powers differential expression analyses for RNA-sequencing and microarray studies." Nucleic acids research 43.7 (2015): e47-e47.



### Disclaimer:
# Hi, my name is Tobias and I hope you can make use
# of this script. Please excuse if the code is not elegant
# as I am a beginner with R. I am thankful for suggestions.
# This script was inspired by the "classical LFQbench papers"
# (Kuharev et al. and Navarro et al.), 
# and serves to acquire improved performance evaluation in benchmarks.
# (Plots for quality assurance and summmary stats for performance description)



### GitHub
# https://github.com/t-jumel/LFQb



### Contact:
# Tobias Jumel
# jumel@mpi-cbg.de
# Shevchenko Lab, MPI-CBG
# 20.03.2023
# Limma contributed by Andre Gohr, MPI-CBG scientific computing
# with advise from Frank Stein, EMBL Heidelberg



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
# Base packages should install automatically, but problems 
# especially with bioconductor packages are not unheard of.



# 2) Select a folder from which input is extracted.
# The folder location is assigned to the variable "folder_input"

# I can recommend "path copy copy" to reformat a windows
# folder location link into an R usable variant.

# This folder must contain exactly one protein group
# and one precursor .tsv matrix with DIA-NN style column names
# with "pg_matrix" and "pr_matrix" in the respective filenames.

# Usage of +MBR requires either the undesired first-pass or 
# second pass matrices to be moved away , e.g. into a subfolder.



# 3) Adjust other variables. 
# This mostly refers to indices of columns containing quantitative values,
# expected log2 fold-change values and filter settings.

# The script is designed to work with 2-4 species in order (human, yeast, ecoli, celegans)
# If your benchmark has different species setup you want to analyse
# you have to search and replace the species term within the code.

# If your samples have just less species but in order
# (human+yeast or human+e.coli) you can execute the script directly,
# with just setting expected log2 fold-change values

# This script assumes human to be stable, yeast up, and both E.coli
# and C.elegans to be downregulated in experimental vs. control condition. 

# In case your samples have exactly opposite up and downregulation trends
# you swap column indices in variables.

# Any other trend deviation requires changes to the code
# deriving the number of true positives etc. from differential expression analysis.



# 4) Execute and find the output folder with csv files and plots
# within the specified input folder.

# For a guide on interpretation see GitHub
# https://github.com/t-jumel/LFQb




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



# Install and load Bioconductor packages automaltically
# Update or install often prohibited by RStudio not allowed to write in folders.
# Try open RStudio "run as administrator" and use .libPaths() to find and manually
# delete folders via windows file explorer to allow fresh re-install. 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)

if (!requireNamespace("limma", quietly = TRUE))
  BiocManager::install("limma", update = FALSE)
library(limma)

if (!requireNamespace("Biobase", quietly = TRUE))
  BiocManager::install("Biobase", update = FALSE)
library(Biobase)




# Variables ---------------------------------------------------------------

# Explanation below at section end for better spacing and overview.


# Set working directory (just pre-typed, might not be used).
# setwd()

# Script version and name used in output filenames.
script_name <- "LFQb"
script_version <- "v3.4"



# Recommend use of "path copy copy".
# See "To run this Script" section for input requirements.
folder_input <- "//fileserver/project_jumel/All_LFQ_benchmarks/External/2021_PXD028735_VanPuyvelde/timsTOFPro_diaPASEF/01_089"



# Naming the 2 different benchmark sample mixtures.
# This is just for record-keeping and preventing manual errors 
# when assigning columns to input categories below.
cond_exp <- "LFQ_A"
cond_ctr <- "LFQ_B"


# From input tables, provide indices of columns
# with quantitative values for the 2 conditions named above.
# Pretyped variants for DIA-NN v1.8 Output layout.
# "Prec" stands for precursor, "Prot" for protein group

# # Pretyped 3 replicates per sample type.
col_exp_Prot <- 6:8
col_ctr_Prot <- 9:11
col_exp_Prec <- 11:13
col_ctr_Prec <- 14:16

# # Pretyped 3 replicates from set of 4 replicates.
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

rep_max <- max(c(length(col_ctr_Prot), length(col_exp_Prot)))
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
plot_res <- 600




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



# Export used variables as separate .csv file for record keeping.
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
# Removal of filetype suffix is performed multiple times
# according to common occurrences.
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
  a.df <- a.df[!grepl("TREMBL", a.df[, "Protein.Names"]),]
  a.df <- a.df[!grepl("SWISS-PROT", a.df[, "Protein.Names"]),]
  a.df <- a.df[!grepl("REFSEQ", a.df[, "Protein.Names"]),]
  a.df <- a.df[!grepl("ENSEMBL", a.df[, "Protein.Names"]),]
  a.df <- a.df[!grepl("H-INV:HIT", a.df[, "Protein.Names"]),]
  a.df <- a.df[!grepl("Streptavidin", a.df[, "Protein.Names"]),]
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






f.species <- function(a.df) {
  
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
  
  
  # Add column listing expected log2 fold-change row-wise
  # to simplify calculations of summary_stats around the log2F fold-change.
    a.df[, "expected_log2FC"] <- NA
  
    a.df[a.df[, "Species"]== "Human", "expected_log2FC"]      <- expFC_human
    a.df[a.df[, "Species"]== "Yeast", "expected_log2FC"]      <- expFC_yeast
    a.df[a.df[, "Species"]== "E.coli", "expected_log2FC"]     <- expFC_ecoli
    a.df[a.df[, "Species"]== "C.elegans", "expected_log2FC"]  <- expFC_celegans
  
  
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
    density <- density(na.omit(data))
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
    density <- density(na.omit(data))
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







# Import, Processing, Filtering of Protein Groups ------------------------------------------------------------------


# Import protein group and precursor matrix. from a selected folder.
Prot0 <- f.import.matrix(folder_input, "pg_matrix", "\t")
Prec0 <- f.import.matrix(folder_input, "pr_matrix", "\t")



## Protein group matrix sanitation and calculations



# For columns with selected quantitative values in DIA-NN matrices,
# reduce text from filepath to replicate name (remove file suffix and filepath)
Prot1 <- Prot0
Prot1 <- f.strip.sample.column.name(Prot1,c(col_ctr_Prot, col_exp_Prot))



# Remove maxQuant contaminant entries, count the removed for record keeping.
Prot1 <- f.remove.MQ.cont(Prot1)
Count_MQcont_Prot <- nrow(Prot0) - nrow(Prot1)



# Replace 0 with NA. "0" is sometimes used as mock value.
# The mock value should not enter any calculation
# and the event should be seen as missing value.(NA)
Prot1[col_ctr_Prot][Prot1[col_ctr_Prot] == 0] <- NA
Prot1[col_exp_Prot][Prot1[col_exp_Prot] == 0] <- NA



# Row-wise summary stats for quant. values within each of the two conditions 
# (experimental vs control)
# Count, Mean, Standard Deviation and CV (Coefficient of Variation)
# These get either "exp" or "ctr" suffix to differentiate the 2 conditions further.
# Add the log2 fold-change
Prot1 <- f.row.stats(Prot1, col_ctr_Prot, col_exp_Prot)
Prot1[, "log2FC"] <- log2(Prot1[, "exp_mean"] / Prot1[, "ctr_mean"])



# Assign Species to rows.
# Remove entries with mixed species.
Prot1 <- f.species(Prot1)
Count_MixSpecies_Prot <- (nrow(Prot0) - Count_MQcont_Prot) - nrow(Prot1)



# Retain entries matching defined data missingness or completeness
# over replicates, typically 50% in both conditions simultaneously.
# Passing entries count as "IDs".
Prot2 <- subset(Prot1, ctr_count >= rep_min & exp_count >= rep_min)



# Continue filtering by precision (CV)
# Passing entries count as "Quant." (quantifications)
Prot3 <- subset(Prot2, ctr_CV <  limit_CV & exp_CV <  limit_CV)




# Differential Expression Analysis by Limma ----------------------------------------

# Limma to get BH adjusted p-values (p_adj) 


# Matrix of log2 transformed quantitative values with proteinn groups as rownames.
assay = as.matrix(Prot3[, c(col_exp_Prot, col_ctr_Prot)]) 
rownames(assay) = Prot3[, 'Protein.Group']
assay = log2(assay)


##List of conditions for all replicates
grps = c(rep('exp', length(col_exp_Prot)), rep('ctr', length(col_ctr_Prot))) 
limma.cond = factor(grps)


# Comparisons to be performed (experimental vs. control)
comparisons <- c('exp - ctr') 


# Create large expression set
limma_data = ExpressionSet(
  assay,
  phenoData = AnnotatedDataFrame(data.frame(CONDITION = grps, 
                                  row.names = colnames(assay))),
  
  featureData = AnnotatedDataFrame(data.frame(PROTEIN_GROUP = Prot3[, 'Protein.Group'],
                                                row.names = rownames(assay)))     
                            )

# Contrast Matrix
contrast.matrix <- model.matrix(~  0 + limma.cond)
colnames(contrast.matrix) <- gsub("limma.cond", "", colnames(contrast.matrix))


#Limma model fit
limma.object <- eBayes(contrasts.fit(
  lmFit(limma_data, design = contrast.matrix),
  makeContrasts(contrasts = comparisons, levels = contrast.matrix)),
trend = TRUE,
robust = TRUE)


# Limma Result Matrix
res_limma = limma::topTable(
  limma.object,
  coef = comparisons[1],
  number = Inf,
  sort.by = 'none',
  adjust.method = "BH",
  confint = TRUE
)


# Add p_adj (limma result) to the main protein group data frame.
Prot3[, "p_adj"] <- res_limma[, 'adj.P.Val']
Prot3[, "p_adj"] <- as.numeric(format(Prot3[, "p_adj"],
                                     scientific = FALSE,
                                     justified = "none"))



# Cleanup environment
rm(assay)
rm(contrast.matrix)
rm(limma_data)
rm(limma.object)
rm(comparisons)
rm(grps)
rm(limma.cond)
rm(limma.rep)
rm(reps)
# rm(res_limma)




# Classification of entries as either
# true positive, true negative, false positive or false negative
# based on the adjusted p-value from limma, the measured log2 fold-change,
# and the expected trend.

# The expectation of human-stable, yeast-upregulated,
# end ecoli, celegans downregulated is hard-coded
# and might have to be adapted for sample layouts.


Prot3[(Prot3[, "Species"] == "Human"
      & (Prot3[, "p_adj"] >= alpha_limma
         | abs(Prot3[, "log2FC"]) <= limit_FC)), "DE_result"] <-
  "true negative"

Prot3[(Prot3[, "Species"] == "Human"
      & Prot3[, "p_adj"] < alpha_limma
      & abs(Prot3[, "log2FC"]) > limit_FC), "DE_result"] <-
  "false positive"



Prot3[(Prot3[, "Species"] == "Yeast"
      & Prot3[, "log2FC"] > +limit_FC
      & Prot3[, "p_adj"] < alpha_limma), "DE_result"] <-
  "true positive"

Prot3[(Prot3[, "Species"] == "Yeast"
      & Prot3[, "log2FC"] < -limit_FC
      & Prot3[, "p_adj"] < alpha_limma), "DE_result"] <-
  "false positive"

Prot3[(Prot3[, "Species"] == "Yeast"
      & (Prot3[, "p_adj"] >= alpha_limma
         | abs(Prot3[, "log2FC"]) <= limit_FC)), "DE_result"] <-
  "false negative"



Prot3[(Prot3[, "Species"] == "E.coli"
      & Prot3[, "log2FC"] < -limit_FC
      & Prot3[, "p_adj"] < alpha_limma), "DE_result"] <-
  "true positive"

Prot3[(Prot3[, "Species"] == "E.coli"
      & Prot3[, "log2FC"] > +limit_FC
      & Prot3[, "p_adj"] < alpha_limma), "DE_result"] <-
  "false positive"

Prot3[(Prot3[, "Species"] == "E.coli"
      & (Prot3[, "p_adj"] >= alpha_limma
         | abs(Prot3[, "log2FC"]) <= limit_FC)), "DE_result"] <-
  "false negative"



Prot3[(Prot3[, "Species"] == "C.elegans"
      & Prot3[, "log2FC"] < -limit_FC
      & Prot3[, "p_adj"] < alpha_limma), "DE_result"] <-
  "true positive"

Prot3[(Prot3[, "Species"] == "C.elegans"
      & Prot3[, "log2FC"] > +limit_FC
      & Prot3[, "p_adj"] < alpha_limma), "DE_result"] <-
  "false positive"

Prot3[(Prot3[, "Species"] == "C.elegans"
      & (Prot3[, "p_adj"] >= alpha_limma
         | abs(Prot3[, "log2FC"]) <= limit_FC)), "DE_result"] <-
  "false negative"







# Processing Precursor Data -----------------------------------------------

# As above for protein groups.
# but without differential expression analysis.

# For columns with qunatitative values in DIA-NN matrices,
# reduce text from filepath to replicate name (remove file suffix and filepath)
Prec1 <- Prec0
Prec1 <- f.strip.sample.column.name(Prec1,c(col_ctr_Prec, col_exp_Prec))



# Remove maxQuant contaminant entries, count the removed for record keeping.
Prec1 <- f.remove.MQ.cont(Prec1)
Count_MQcont_Prec <- nrow(Prec0) - nrow(Prec1)



# Replace 0 with NA. "0" is sometimes used as mock value.
# The mock value should not enter any calculation
# and the event should be seen as missing value.(NA)
Prec1[col_ctr_Prec][Prec1[col_ctr_Prec] == 0] <- NA
Prec1[col_exp_Prec][Prec1[col_exp_Prec] == 0] <- NA



# Row-wise summary stats for quant. values within each of the two conditions 
# (experimental vs control)
# Count, Mean, Standard Deviation and CV (Coefficient of Variation)
# These get either "exp" or "ctr" suffix to differentiate the 2 conditions further.
# Add the log2 fold-change
Prec1 <- f.row.stats(Prec1, col_ctr_Prec, col_exp_Prec)
Prec1[, "log2FC"] <- log2(Prec1[, "exp_mean"] / Prec1[, "ctr_mean"])



# Assign Species to rows.
# Remove entries with mixed species.
Prec1 <- f.species(Prec1)
Count_MixSpecies_Prec <- (nrow(Prec0) - Count_MQcont_Prec) - nrow(Prec1)



# Treat Precursors as unfiltered by data completeness and CV,
# but generate an alternate filtered version analogue to Prot2 and Prot3.
Prec2 <- subset(Prec1, ctr_count >= rep_min & exp_count >= rep_min)
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
  skewness(na.omit(Prec3[which(Prec3[, "Species"] == "E.coli"), "log2FC"])),
  skewness(na.omit(Prec3[which(Prec3[, "Species"] == "Human"), "log2FC"])),
  skewness(na.omit(Prec3[which(Prec3[, "Species"] == "Yeast"), "log2FC"])),
  skewness(na.omit(Prec3[which(Prec3[, "Species"] == "C.elegans"), "log2FC"])),
  skewness(na.omit(Prot3[which(Prot3[, "Species"] == "E.coli"), "log2FC"])),
  skewness(na.omit(Prot3[which(Prot3[, "Species"] == "Human"), "log2FC"])),
  skewness(na.omit(Prot3[which(Prot3[, "Species"] == "Yeast"), "log2FC"])),
  skewness(na.omit(Prot3[which(Prot3[, "Species"] == "C.elegans"), "log2FC"])))

asymmetry[, "Kurtosis"] <- c(
  kurtosis(na.omit(Prec3[which(Prec3[, "Species"] == "E.coli"), "log2FC"])),
  kurtosis(na.omit(Prec3[which(Prec3[, "Species"] == "Human"), "log2FC"])),
  kurtosis(na.omit(Prec3[which(Prec3[, "Species"] == "Yeast"), "log2FC"])),
  kurtosis(na.omit(Prec3[which(Prec3[, "Species"] == "C.elegans"), "log2FC"])),
  kurtosis(na.omit(Prot3[which(Prot3[, "Species"] == "E.coli"), "log2FC"])),
  kurtosis(na.omit(Prot3[which(Prot3[, "Species"] == "Human"), "log2FC"])),
  kurtosis(na.omit(Prot3[which(Prot3[, "Species"] == "Yeast"), "log2FC"])),
  kurtosis(na.omit(Prot3[which(Prot3[, "Species"] == "C.elegans"), "log2FC"])))


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








# Counts of protein groups and precursors

Counts <- data.frame(matrix(ncol = 0, nrow = 5))
rownames(Counts) <- c("Total entries",
                      "Removed MaxQuant contaminants",
                      "Removed mixed species entries",
                      "Removed for data incompletness",
                      "Removed for quantitative imprecision"
                      )


Counts[, "Precursors"] <- c(nrow(Prec0),
                              Count_MQcont_Prec,
                              Count_MixSpecies_Prec,
                              nrow(Prec1) - nrow(Prec2),
                              nrow(Prec2) - nrow(Prec3))

Counts[, "Protein Groups"] <- c(nrow(Prot0),
                           Count_MQcont_Prot,
                           Count_MixSpecies_Prot,
                           nrow(Prot1) - nrow(Prot2),
                           nrow(Prot2) - nrow(Prot3))


write.csv(Counts,
          row.names = TRUE,
          paste0(folder_output, "/", "counts_",
                 basename(dirname(folder_output)),"_",
                 script_version,"_",
                 rep_min,"of", rep_max,
                 "_CV", limit_CV,
                 "_FC",limit_FC,
                 "_a",alpha_limma, ".csv"))





# Main summary stats.

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
          100 * length(which(Prot3[, "DE_result"] == "false positive")) /
            (length(which(
              Prot3[, "DE_result"] == "false positive"
            ))
            + length(which(
              Prot3[, "DE_result"] == "true positive"
            ))))
  
  
  # true positive count
  summary_stats[, "TP"] <-
    length(which(Prot3[, "DE_result"] == "true positive"))
  
  
  
  
  summary_stats[, "Sensitivity"] <-
    round(digits = 2, 100 * length(which(Prot3[, "DE_result"] == "true positive")) /
            (length(which(
              Prot3[, "DE_result"] == "true positive"
            )) + length(which(
              Prot3[, "DE_result"] == "false negative"
            ))))
  
  summary_stats[, "Specificity"] <-
    round(digits = 2, 100 * length(which(Prot3[, "DE_result"] == "true negative")) /
            (length(which(
              Prot3[, "DE_result"] == "true negative"
            )) + length(which(
              Prot3[, "DE_result"] == "false positive"
            ))))
  
  
  
  
  
  # Number of confidently detected and quantified protein groups,
  # after filtering for missingness and missingness + CV.
  summary_stats[, "Prot_ID"] <- nrow(Prot2)
  summary_stats[, "Prot_Quant"] <- nrow(Prot3)
  
  
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
    round(digits = 2,  mean(c(Prot3[, "ctr_CV"], Prot3[, "exp_CV"])))
  summary_stats[, "Prot_CV_Median"] <-
    round(digits = 2, median(c(Prot3[, "ctr_CV"], Prot3[, "exp_CV"])))
  
  
  # Accuracy is average distance of measured to expected log2FC.
  summary_stats[, "Prot_Accuracy"] <-
    round(digits = 2,
          mean(abs(Prot3[, "log2FC"] - Prot3[, "expected_log2FC"])))
  
  
  # Trueness is distance between measurement medians 
  # and respective expected fold-changes.
  # (cumulative as shifts are important to spot)
  # Several conditions can lead to a high value
  # e.g. ratio compression or erroneous normalization or even both together.
  summary_stats[, "Prot_Trueness"] <-
    round(
      digits = 2,
      sum(na.rm = TRUE,
          c(abs(median(Prot3[Prot3[, "Species"] == "Human"  , "log2FC"]) - expFC_human),
            abs(median(Prot3[Prot3[, "Species"] == "Yeast"  , "log2FC"]) - expFC_yeast),
            abs(median(Prot3[Prot3[, "Species"] == "E.coli" , "log2FC"]) - expFC_ecoli),
            abs(median(Prot3[Prot3[, "Species"] == "C.elegans"  , "log2FC"]) - expFC_celegans)
          )))
  
  
  # Dispersion is the average distance of individual measurements around the respective median.
  summary_stats[, "Prot_Dispersion"] <-
    round(digits = 2,
          mean(na.rm = TRUE,
               c(
                 abs(Prot3[Prot3[, "Species"] == "Human"   , "log2FC"] - median(Prot3[Prot3[, "Species"] == "Human"  , "log2FC"])),
                 abs(Prot3[Prot3[, "Species"] == "Yeast"   , "log2FC"] - median(Prot3[Prot3[, "Species"] == "Yeast"  , "log2FC"])),
                 abs(Prot3[Prot3[, "Species"] == "E.coli"  , "log2FC"] - median(Prot3[Prot3[, "Species"] == "E.coli" , "log2FC"])),
                 abs(Prot3[Prot3[, "Species"] == "C.elegans"  , "log2FC"] - median(Prot3[Prot3[, "Species"] == "C.elegans" , "log2FC"]))
               )))
  
  
  # Other than TP, continued confusion matrix of 
  # differential expression result interpretation.
  summary_stats[, "FP"] <-
    length(which(Prot3[, "DE_result"] == "false positive"))
  summary_stats[, "TN"] <-
    length(which(Prot3[, "DE_result"] == "true negative"))
  summary_stats[, "FN"] <-
    length(which(Prot3[, "DE_result"] == "false negative"))
  
  
  
  
  
  
  
  
  # Stats partially repeated on precursor level as for protein groups above.
  summary_stats[, "Prec_ID"] <- nrow(Prec2)
  summary_stats[, "Prec_Quant"] <- nrow(Prec3)
  
  
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
                                              mean(c(Prec3[, "ctr_CV"], Prec3[, "exp_CV"])))
  summary_stats[, "Prec_CV_Median"] <-  round(digits = 2,
                                              median(c(Prec3[, "ctr_CV"], Prec3[, "exp_CV"])))
  
  
  # Accuracy is average distance of measured to expected log2FC.
  summary_stats[, "Prec_Accuracy"] <-
    round(digits = 2,
          mean(abs(Prec3[, "log2FC"] - Prec3[, "expected_log2FC"])))
  
  
  # Trueness is distance between measurement medians 
  # and respective expected fold-changes (cumulative as shifts are important to spot)
  summary_stats[, "Prec_Trueness"] <-
    round(
      digits = 2,
      sum(na.rm = TRUE,
          c(abs(median(Prec3[Prec3[, "Species"] == "Human"  , "log2FC"]) - expFC_human),
            abs(median(Prec3[Prec3[, "Species"] == "Yeast"  , "log2FC"]) - expFC_yeast),
            abs(median(Prec3[Prec3[, "Species"] == "E.coli"  , "log2FC"]) - expFC_ecoli),
            abs(median(Prec3[Prec3[, "Species"] == "C.elegans"  , "log2FC"]) - expFC_celegans)
          )))
  
  
  # Dispersion is the average distance of individual measurements around the respective median.
  summary_stats[, "Prec_Dispersion"] <-
    round(digits = 2,
          mean(na.rm = TRUE,
               c(
                 abs(Prec3[Prec3[, "Species"] == "Human"  , "log2FC"] - median(Prec3[Prec3[, "Species"] == "Human"  , "log2FC"])),
                 abs(Prec3[Prec3[, "Species"] == "Yeast"  , "log2FC"] - median(Prec3[Prec3[, "Species"] == "Yeast"  , "log2FC"])),
                 abs(Prec3[Prec3[, "Species"] == "E.coli"  , "log2FC"] - median(Prec3[Prec3[, "Species"] == "E.coli"  , "log2FC"])),
                 abs(Prec3[Prec3[, "Species"] == "C.elegans"  , "log2FC"] - median(Prec3[Prec3[, "Species"] == "C.elegans"  , "log2FC"]))
               )))
  
  
  
  # Folder that contained the analysed pg and pc matrix
  summary_stats[, "folder_input"] <- folder_input

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
  " \n Sens. ",
  summary_stats[, "Sensitivity"],
  " %, ",
  "Spec. ",
  summary_stats[, "Specificity"],
   " %\n",
  summary_stats[, "TP"],
  " true positives, ",
  summary_stats[, "deFDR"],
  "% deFDR")

subtitle_Prot_density <- paste0(
  summary_stats[, "Prot_Quant"],
  " Quant. (",
  round(nrow(Prot3) / nrow(Prot2) * 100, digits = 0),
  "%), ",
  summary_stats[, "Prot_ID"],
  " IDs",
  " \n Asym. E.coli ",
  summary_stats[, "Prot_Asymmetry_E.coli"],
  ", ",
  "Asym. Yeast ",
  summary_stats[, "Prot_Asymmetry_Yeast"],
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
    " \n Asym. E.coli ",
    summary_stats[, "Prec_Asymmetry_E.coli"],
    ", ",
    "Asym. Yeast ",
    summary_stats[, "Prec_Asymmetry_Yeast"],
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


# Plot_CVbox ------------------------------------------------------------


f.p.CV_box <- function(a.df,
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
    
    geom_boxplot(
      key_glyph = "blank",
      outlier.colour = NA,
      outlier.shape = 16,
      outlier.size = 1,
      outlier.alpha = 0.3,
      alpha = 0.7,
      notch = FALSE,
      fatten = 1,
      show.legend = FALSE
    )+
    
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
    
    coord_cartesian(ylim = c(0, 20)) +
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
plot_Prot_CV_box <- f.p.CV_box(Prot3, subtitle_Prot, caption_Prot)
plot_Prec_CV_box <- f.p.CV_box(Prec3, subtitle_Prec, caption_Prec)


# Export plots
f.export.plot(plot_Prot_CV_box, "Prot_CV_box", plot_width, plot_height, plot_res)
f.export.plot(plot_Prec_CV_box, "Prec_CV_box", plot_width, plot_height, plot_res)






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
                                 subtitle_Prot_density,
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



# Plot_FCbox ------------------------------------------------------------


f.p.FC_box <- function(a.df,
                       a.subtitle,
                       a.caption) {
  ggplot((
    
      a.df
  ),
  aes(x = Species,
      y = abs(a.df[, "log2FC"]-a.df[, "expected_log2FC"]),
      fill = Species)) +
    
    geom_boxplot(
      key_glyph = "blank",
      outlier.colour = NA,
      outlier.shape = 16,
      outlier.size = 1,
      outlier.alpha = 0.3,
      alpha = 0.7,
      notch = FALSE,
      fatten = 1,
      show.legend = FALSE
    )+
    
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
    
    coord_cartesian(
      ylim = c(0, 1.0)
      ) +
    scale_y_continuous(
      name = "|expected - measured log2 fold-change|",
      #breaks = c(0, 5, 10, 15, 20, 25, 30),
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
plot_Prot_FC_box <- f.p.FC_box(Prot3, subtitle_Prot, caption_Prot)
plot_Prec_FC_box <- f.p.FC_box(Prec3, subtitle_Prec, caption_Prec)


# Export plots
f.export.plot(plot_Prot_FC_box, "Prot_FC_box", plot_width, plot_height, plot_res)
f.export.plot(plot_Prec_FC_box, "Prec_FC_box", plot_width, plot_height, plot_res)









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