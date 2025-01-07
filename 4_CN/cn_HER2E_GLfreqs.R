# Script: calculate frequencies of gain/loss alterations per probe
# Purpose: 
# 1. calculate frequencies of gain/loss alterations per probe

# input file from cn_HER2Eprocessing.R

# TODO:
# - 

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run (for plot naming)
cohort <- "SCANB"

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/4_CN/processed/",sep="")
dir.create(data.path)

# set/create output directory for plots
output.path <- "output/plots/4_CN/"
dir.create(output.path)

#packages
source("scripts/4_CN/src/cn_functions.R")
library(ggplot2)
library(tidyverse)

################################################################################
# loading data
################################################################################

# this file has NA probes filtered out -> correct? only if i plot genes and not probes otherwise they can stay in
gainloss.cn.scanb.freqs <- loadRData("data/SCANB/4_CN/processed/CN_gainloss_genpos_genmap.RData") %>% 
  dplyr::select(-c(Gene_symbol))

################################################################################
# calc. the gain/loss freq. for HER2E 
################################################################################

# calc loss/gain freqs per group
gainloss.cn.scanb.freqs$freqloss.HER2E <- apply(
  gainloss.cn.scanb.freqs[,5:ncol(gainloss.cn.scanb.freqs)], 1, function(x) (
    length(which(x<=-1))/ncol(gainloss.cn.scanb.freqs[,5:ncol(gainloss.cn.scanb.freqs)]))*-100) # i add a minus to make it easier for plotting

gainloss.cn.scanb.freqs$freqgain.HER2E <- apply(
  gainloss.cn.scanb.freqs[,5:ncol(gainloss.cn.scanb.freqs)], 1, function(x) (
    length(which(x>=1))/ncol(gainloss.cn.scanb.freqs[,5:ncol(gainloss.cn.scanb.freqs)]))*100)

#head(gainloss.cn.scanb.freqs)

# remove the indiv. sample data
gainloss.cn.scanb.freqs <- gainloss.cn.scanb.freqs %>% 
  relocate(freqgain.HER2E, .after=Genome_pos) %>% 
  relocate(freqloss.HER2E, .after=Genome_pos) 
gainloss.cn.scanb.freqs[,7:ncol(gainloss.cn.scanb.freqs)] <- NULL

#head(gainloss.cn.scanb.freqs)

# save
save(gainloss.cn.scanb.freqs,
     file = paste(data.path,"CN_gainloss_frequencies.RData",sep=""))
