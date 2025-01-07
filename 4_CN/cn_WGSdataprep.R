# Script: Preparing analyses based on WGS data
# includes generating cn summary matrices 

#TODO: 
# plot to test if the probe states and segment states match

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

# plot
plot.list <- list() # object to store plots; note: if the output is not in string format use capture.output()
plot.file <- paste(output.path,cohort,"_HER2n_WGS.pdf",sep = "")

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

#View(ascat.files[[1]])

# select only the segments df and remove sample column
segment.files <- lapply(ascat.files, function(x) { 
  x <- x[["segments"]] %>% dplyr::select(-c(sample))
  x$nTotal <- x$nMinor + x$nMajor
  return(x) # add nTotal column
  })

#str(segment.files) # not named yet
#View(segment.files[1])
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

save(ploidy.df,file = "./data/SCANB/4_CN/processed/HER2Esample_ploidy.RData")

# probe anno
probe.df <- read.delim(
  "./data/SCANB/3_genomic/raw/CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+/ascat/SnpGcCorrections.tsv",sep="\t")[1:3] %>% 
  mutate(Chr = gsub("chr","",Chr)) %>% 
  dplyr::rename("ProbeID"=X)
#str(probe.df)

#######################################################################
# functions (move to CN functions later)
#######################################################################
#i=1
# function that that creates a matrix (probeID | Chr | Position | CNstate_sample) for all samples
make_finmatrix <- function(probe.df, segment.files, ploidy.df, matrix.type) {
  
  fin.list <- list()
  
  # for each file (sample)
  for(i in c(1:length(segment.files))) { 
    
    # sample data
    sample.ID <- names(segment.files[i]) 
    sample.data <- segment.files[[sample.ID]] 
    sample.ploidy <- unlist(ploidy.df[ploidy.df$Sample==sample.ID,]$ploidy)
    
    # create a matrix (probeID | Chr | Position | CNstate_sample) for a single sample
    sample.matrix <- make_sampmatrix(probe.df, sample.data, sample.ID, matrix.type, sample.ploidy)
    
    # list with sample matrices
    fin.list <- append(fin.list, list(sample.matrix))
  }

  # merge fin.list into final matrix
  fin.matrix <- fin.list %>% purrr::reduce(left_join,by=c("ProbeID","Chr","Position"))
  
  return(fin.matrix)
}

#######################################################################

## function that that creates a matrix (probeID | Chr | Position | CNstate_sample) for a single sample
make_sampmatrix <- function(probe.df, sample.data, sample.ID, matrix.type, sample.ploidy) {
  
  chr.set <- unique(sample.data$chr)
  res.list <- list()
  
  for(j in c(1:length(chr.set))) { 
    # chromosome
    chr <- chr.set[j] #j
  # relevant probes
    chr.probes <- probe.df[probe.df$Chr == chr,]
    # relevant segments
    chr.segments <- sample.data[sample.data$chr == chr,]
    # get chr matrix
    chr.matrix <- make_chrmatrix(
      chr.probes, chr.segments, matrix.type, sample.ploidy)
    
    # save in list
    res.list <- append(res.list,list(chr.matrix))
  }
  
  #make one dataframe out of the list
  sample.matrix <- do.call(rbind, res.list)
  
  # rename cn_state to the sample id
  sample.matrix <- sample.matrix %>% dplyr::rename(!!sample.ID := CN_state)
  
  return(sample.matrix)
}

#######################################################################

# function that creates a CN matrix (probeID | Chr | Position | CNstate_sample) for a single chromosome with a method defining what type of matrix is created (total, major, minor)
make_chrmatrix <- function(chr.probes, chr.segments, matrix.type, sample.ploidy) {
  
  # add a column to probe file that is to be filled with the CN state for each probe
  chr.matrix <- chr.probes
  chr.matrix$CN_state <- NA 
  
  # for each chromosome segment
  for(k in c(1:nrow(chr.segments))) {
    
    # get sample info
    seg.data <- chr.segments[k,] 
    seg.end <- seg.data[["endpos"]]
    seg.start <- seg.data[["startpos"]]
    
    # define calculations for each matrix type
    # major CN state
    if (matrix.type == "major") {
        # assign the CN_state to the probes that correspond to that segment
        chr.matrix[which(chr.matrix$Position >= seg.start & 
                           chr.matrix$Position <= seg.end),]$CN_state <- seg.data[["nMajor"]]
    # minor CN state
      } else if (matrix.type == "minor") {
        # assign the CN_state to the probes that correspond to that segment
        chr.matrix[which(chr.matrix$Position >= seg.start & 
                           chr.matrix$Position <= seg.end),]$CN_state <- seg.data[["nMinor"]]
    # total CN state  
      } else if (matrix.type == "total") {
        # assign the CN_state to the probes that correspond to that segment
        chr.matrix[which(chr.matrix$Position >= seg.start & 
                           chr.matrix$Position <= seg.end),]$CN_state <- seg.data[["nTotal"]]
    # gain or loss CN state
      } else if (matrix.type == "gainloss") {
        if (seg.data[["nTotal"]]  >= (sample.ploidy + 0.6)) {
          state.gl <- 1
          } else if (seg.data[["nTotal"]] <= (sample.ploidy - 0.6)) {
            state.gl <- -1
          } else {
            state.gl <- 0
        }
        # assign the CN_state to the probes that correspond to that segment
        chr.matrix[which(chr.matrix$Position >= seg.start & 
                           chr.matrix$Position <= seg.end),]$CN_state <- state.gl
    # amplification CN state
      } else if (matrix.type == "amplification") {
        if (seg.data[["nTotal"]]  >= (sample.ploidy*4)) { # focal amp
          state.amp <- 2
          } else if (seg.data[["nTotal"]] >= (sample.ploidy*2)) { # amp
            state.amp <- 1
          } else {
            state.amp <- 0 # no amp
        }
        # assign the CN_state to the probes that correspond to that segment
        chr.matrix[which(chr.matrix$Position >= seg.start & 
                           chr.matrix$Position <= seg.end),]$CN_state <- state.amp
    # LOH CN state
      } else if (matrix.type == "LOH") {
        if (seg.data[["nMinor"]]  == 0) { # LOH
          state.loh <- 1
          } else {
            state.loh <- 0 # no LOH
          }
        # assign the CN_state to the probes that correspond to that segment
        chr.matrix[which(chr.matrix$Position >= seg.start & 
                           chr.matrix$Position <= seg.end),]$CN_state <- state.loh
      } else { print("Check matrix type input") }
  } 
  return(chr.matrix)
}

################ test

#res <- make_chrmatrix(chr.probes, chr.segments, matrix.type="major", sample.ploidy)
#res2 <- make_sampmatrix(probe.df, sample.data, sample.ID, matrix.type="major", sample.ploidy)
# test for one file
#out <- make_finmatrix(probe.df, segment.files, ploidy.df, matrix.type="amplification")

################ 
vars <- c("major", "minor", "total", "amplification", "gainloss", "LOH")
for (i in c(1:length(vars))) {
  matrix.type <- vars[i]
  # get matrix
  res <- make_finmatrix(probe.df = probe.df, 
                        segment.files = segment.files, 
                        ploidy.df = ploidy.df, 
                        matrix.type = matrix.type)
  res$Chr <- as.numeric(res$Chr)
  # save
  save(res,file=paste(data.path,"CN_",matrix.type,".RData",sep=""))
}
