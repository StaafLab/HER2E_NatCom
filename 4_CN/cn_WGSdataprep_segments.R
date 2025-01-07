# Script: Preparing analyses based on WGS data: assigning CN states to segments

#TODO: 


# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "SCANB" # SCANB 

# set/create output directory for plots
output.path <- "output/plots/4_CN/"
dir.create(output.path)

# set/create output directory for data
data.path <- "data/SCANB/4_CN/raw/"
dir.create(data.path)

#packages
source("scripts/4_CN/src/cn_functions.R")
library(ggplot2)
library(tidyverse)
library(reshape2)
#library(data.table)
library(purrr)
library(readxl)

#######################################################################
#######################################################################

# get normal IDs
id.key <- read_excel("data/SCANB/3_genomic/raw/JVCimpression2022-11-09_ERpos_WGS_Batch1-3.xlsx") %>% 
  dplyr::select(c(SENT.TUMOR,SENT.TUMOR.aliquot,TUMOR.alias)) %>% 
  mutate(SENT.TUMOR.aliquot = gsub("\\.","_",SENT.TUMOR.aliquot))

# compile all files into one R object
# ascat data file paths
temp <- list.files(
  path="./data/SCANB/4_CN/raw/to_lennart/",
  pattern="*.RData", full.names = TRUE, recursive = TRUE)
# load the files
ascat.files <- lapply(temp, loadRData)

# select only the segments df and remove sample column
segment.files <- lapply(ascat.files, function(x) { 
  x <- x[["segments"]] %>% dplyr::select(-c(sample))
  x$nTotal <- x$nMinor + x$nMajor
  return(x) # add nTotal column
  })

# convert the file names to S identifiers
# sample ids used in the ascat files
ascat.ids <- as.data.frame(unlist(lapply(temp, function(x) {
  gsub("_vs.*","",
       gsub("./data/SCANB/4_CN/raw/to_lennart//ascat.","",x))}))) %>% dplyr::rename(Ascat.id=1)
# 1. convert epb IDs to normal sample IDs
ascat.ids$Sample <- id.key$SENT.TUMOR[match(
  ascat.ids$Ascat.id, id.key$SENT.TUMOR.aliquot)] 
# 2. convert the other IDs to normal sample IDs
ascat.ids$Sample2 <- id.key$SENT.TUMOR[match(
  ascat.ids$Ascat.id, id.key$TUMOR.alias)] # 2. convert the other IDs to normal sample IDs
ascat.ids$Sample <- ifelse(is.na(ascat.ids$Sample), ascat.ids$Sample2, ascat.ids$Sample)
ascat.ids$Sample2 <- NULL

# name the segments dataframes
names(segment.files) <- ascat.ids$Sample

# ploidy file
ploidy.list <- lapply(ascat.files, function(x) { 
  ploidy <- x[["ploidy"]]
  sample <- names(ploidy)
  ploidy <- unname(ploidy) 
  return(list(sample,ploidy)) 
}) 
ploidy.df <- as.data.frame(do.call(rbind,ploidy.list)) %>% 
  dplyr::rename(Sample=V1,ploidy=V2)
# name to match the segment sample ids
ploidy.df$Sample <- ascat.ids$Sample[match(
  ploidy.df$Sample, ascat.ids$Ascat.id)] 

#save(ploidy.df,file = "./data/SCANB/4_CN/processed/HER2Esample_ploidy.RData")

#######################################################################
#######################################################################

all.samples.list <- list() # list of sample matrices looking like this: seg, amp, gl, loh

# for each file (sample)
pb = txtProgressBar(min = 0, max = length(segment.files), initial = 0, style = 3)
for(j in c(1:length(segment.files))) { 
  setTxtProgressBar(pb,j)
  
  # sample data
  sample.ID <- names(segment.files[j]) 
  sample.data <- segment.files[[sample.ID]] 
  sample.ploidy <- unlist(ploidy.df[ploidy.df$Sample==sample.ID,]$ploidy)
  
  # result cols
  sample.data$GainLoss <- NA 
  sample.data$Amp <- NA 
  sample.data$LOH <- NA 
  
  for(i in 1:nrow(sample.data)) { # loop through the segments
    
    # define calculations for each matrix type
    # gain or loss CN state
    if (sample.data$nTotal[i]  >= (sample.ploidy + 0.6)) {
      state.gl <- 1 # gain
    } else if (sample.data$nTotal[i] <= (sample.ploidy - 0.6)) {
      state.gl <- -1 # loss
    } else {
      state.gl <- 0 # neither
    }
    # assign the CN_state to the segment
    sample.data$GainLoss[i] <- state.gl 
    
    # amplification CN state
    if (sample.data$nTotal[i]  >= (sample.ploidy*4)) { 
      state.amp <- 2 # focal amp
    } else if (sample.data$nTotal[i] >= (sample.ploidy*2)) { # amp
      state.amp <- 1 # amp
    } else {
      state.amp <- 0 # no amp
    }
    # assign the CN_state to the segment
    sample.data$Amp[i] <- state.amp 
    
    # LOH CN state
    if (sample.data$nMinor[i] == 0) { # LOH
      state.loh <- 1 # LOH
    } else {
      state.loh <- 0 # no LOH
    }
    # assign the CN_state to the probes that correspond to that segment
    sample.data$LOH[i] <- state.loh 
  }
  
  all.samples.list <- append(all.samples.list,list(sample.data))
  close(pb)
}

names(all.samples.list) <- names(segment.files)

save(all.samples.list,file = "./data/SCANB/4_CN/processed/Segment_CN_states.RData")

