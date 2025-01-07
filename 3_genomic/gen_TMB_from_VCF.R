# Script: extract TMB from VCF files

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "COMBINED" 

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/3_genomic/processed/",sep="")
dir.create(data.path)

#packages
source("scripts/3_genomic/src/gen_functions.R")
source("scripts/4_CN/src/cn_functions.R")
library(ggplot2)
library(tidyverse)
library(readxl)
library(GenVisR)
library(reshape2)
library(ggstatsplot)
#library(data.table)
library(grid)
library(gridExtra)
library(R.utils)
library(vcfR)
library(car)

#######################################################################
# load TMB data (extract from VCF files)
#######################################################################

# scanb
# ID key file
id.key <- as.data.frame(read_excel("./data/SCANB/3_genomic/raw/HER2_enriched_June23_ForJohan.xlsx", sheet = "Samples")) %>% 
  dplyr::select(c("Sample","Tumour")) %>% dplyr::rename(sample=Sample)

# get full list from VCF files

# compile all files into one R object
caveman.files <- list.files(
  path="./data/SCANB/3_genomic/raw/vcf/",
  pattern="*.annot.muts.vcf.gz$", full.names = TRUE, recursive = TRUE)
pindel.files <- list.files(
  path="./data/SCANB/3_genomic/raw/vcf/",
  pattern="*.annot.vcf.gz$", full.names = TRUE, recursive = TRUE)

# load the files
caveman.files <- lapply(caveman.files, function(x) {
  # load
  file <- read.vcfR(x)
  sample <- gsub(".*caveman/(.+)_vs.*", "\\1", x)
  file.dat <- as.data.frame(file@fix) %>% 
    # make CLPM and ASMD into own column
    mutate(ASMD = as.numeric(gsub(".*ASMD=([0-9.]+);.*", "\\1", INFO))) %>% 
    mutate(CLPM = as.numeric(gsub(".*CLPM=([0-9.]+);.*", "\\1", INFO))) %>% 
    filter(FILTER=="PASS" & CLPM==0 & ASMD >=140) #Caveman counts (CLPM=0,ASMD>=140)_final
  return(c(sample, nrow(file.dat)))
})  
cv <- as.data.frame(do.call(rbind,caveman.files))
names(cv) <- c("sample","caveman_count")

pindel.files <- lapply(pindel.files, function(x) {
  # load
  file <- read.vcfR(x)
  sample <- gsub(".*pindel/(.+)_vs.*", "\\1", x)
  file.dat <- as.data.frame(file@fix) %>% 
    # make repeats column
    mutate(Repeats = as.numeric(gsub(".*REP=([0-9]+);.*", "\\1", INFO))) %>% 
    filter(FILTER=="PASS" & QUAL>=250 & Repeats<10) #Pindel counts (QUAL>=250,Repeats<10)_final
  return(c(sample, nrow(file.dat)))
})

pd <- as.data.frame(do.call(rbind,pindel.files))
names(pd) <- c("sample","pindel_count")

# merge and correct IDs 
res <- merge(cv,pd,by="sample")

# correct ids now
res$sampleID <- id.key$sample[match(res$sample,id.key$Tumour)]
res$sample <- NULL
res$PAM50 <- rep("HER2E", length(sample))
res$N_mut <- as.numeric(res$caveman_count) + as.numeric(res$pindel_count)
scanb.muts <- res

#save(scanb.muts, file="./data/SCANB/3_genomic/processed/mutation_counts.RData")