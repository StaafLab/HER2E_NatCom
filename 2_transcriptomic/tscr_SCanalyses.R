# Script: Danenberg immune part in metabric continued
# Author: Lennart Hohmann
# Date: 01.06.2024
#-------------------
# empty environment 
rm(list=ls())
# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")
# cohort
cohort <- "METABRIC"
#-------------------
# packages
source("./scripts/src/general_functions.R")
source("./scripts/2_transcriptomic/src/tscr_functions.R")
if (!require("pacman")) install.packages("pacman")
#pacman::p_load()
#-------------------
# set/create output directories
# for plots
output.path <- "./output/plots/2_transcriptomic/"
dir.create(output.path)
# for data
data.path <- "./data/METABRIC/2_transcriptomic/processed/"
dir.create(data.path)
#-------------------
# input paths
infile.1 <- "./data/METABRIC/1_clinical/raw/Merged_annotations.RData"
infile.2 <- "./data/METABRIC/2_transcriptomic/processed/SingleCells_PTcounts.RData"
infile.3 <- "./data/METABRIC/2_transcriptomic/processed/SingleCells_PTprop.RData"
infile.4 <- "./data/METABRIC/2_transcriptomic/processed/SingleCells_PTmedians.RData"
infile.5 <- "./"

# output paths
plot.file <- paste0(output.path,cohort,"_singlecellimmune.pdf")
txt.file <- paste0(output.path,cohort,"_singlecellimmune.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

anno <- loadRData(infile.1)
counts.dat <- loadRData(infile.2)
prop.dat <- loadRData(infile.3)

#######################################################################
# pt cell proportions against mRNA expression
#######################################################################



