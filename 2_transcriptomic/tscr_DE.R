# Script: differential expression analysis; cohorts=SCANB, METABRIC
#TODO: 
# Filter based on diff in expression absolute (add column)

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "METABRIC" # SCANB or METABRIC

# set/create output directory for plots
output.path <- "output/plots/2_transcriptomic/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/2_transcriptomic/processed/",sep="")
dir.create(data.path)

#packages
source("scripts/2_transcriptomic/src/tscr_functions.R")
library(ggplot2)
library(tidyverse)
library(matrixStats)
library(pheatmap)
#library(Hmisc)
library(VennDiagram)
library(readxl)
library(ggfortify)
library(janitor)
library(biomaRt)

#######################################################################
# 2. Cohort-specific data preprocessing including selection of  
# the clinical ER+Her2- subtyped samples
#######################################################################

# for SCANB
if (cohort=="SCANB") {
    
    # load annotation data and select subgroup data
    anno <- loadRData(file="./data/SCANB/1_clinical/processed/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data_ERpHER2n.RData") %>%
      dplyr::rename(sampleID = GEX.assay, PAM50 = NCN.PAM50)
    
    # load gex data
    gex.data <- scanb_gex_load(gex.path = "data/SCANB/2_transcriptomic/raw/genematrix_noNeg.Rdata", geneanno.path = "data/SCANB/1_clinical/raw/Gene.ID.ann.Rdata", ID.type = "Gene.Name") %>% 
        dplyr::select(any_of(anno$sampleID)) %>% # select subgroup gex 
        select_if(~ !any(is.na(.))) # otherwise error when scaling

    # log transform FPKM data
    gex.data <- as.data.frame(log2(gex.data + 1))
    # z-transform
    gex.data <- as.data.frame(t(apply(gex.data, 1, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y))))) # for some rows there may be 0 variance so i have to handle these cases
    
#---------------------------------------------------------------------#
    
} else if (cohort=="METABRIC") {
    
    # load annotation data
    load("data/METABRIC/1_clinical/raw/Merged_annotations.RData")
    
    # extract relevant variables
    anno <- loadRData("data/METABRIC/1_clinical/processed/Merged_annotations_ERpHER2n.RData") %>% 
      dplyr::rename(sampleID=METABRIC_ID) # rename to match SCANB variables
    
    # load and select subgroup data
    gex.data <- metabric_gex_load("./data/METABRIC/2_transcriptomic/raw/data_mRNA_median_all_sample_Zscores.txt",ID.type = "Hugo_Symbol") %>% 
        dplyr::select(any_of(anno$sampleID)) %>% 
        mutate_all(function(x) as.numeric(x))
    
    # exclude samples from anno without associated gex data
    anno <- anno %>% 
        filter(sampleID %in% colnames(gex.data))
}

#######################################################################
# 4. Perform t-test for each gene
#######################################################################

# sampleIDs per pam50 subtype
Her2.cols <- anno %>% filter(PAM50=="Her2") %>% pull(sampleID)
LumA.cols <- anno %>% filter(PAM50=="LumA") %>% pull(sampleID)
LumB.cols <- anno %>% filter(PAM50=="LumB") %>% pull(sampleID)

# initialize vector of stored p-values and expression differences (based on mean comparison)
Her2.LumA.pval <- rep(NA,nrow(gex.data))
Her2.LumB.pval <- rep(NA,nrow(gex.data))
Her2.LumA.diff <- rep(NA,nrow(gex.data))
Her2.LumB.diff <- rep(NA,nrow(gex.data))

# loop through the genes
pb = txtProgressBar(min = 0, max = nrow(gex.data), initial = 0, style = 3)
for (i in 1:nrow(gex.data)) { #nrow(gex.data)
    setTxtProgressBar(pb,i)
    
    # set vars
    her2.data <- as.numeric(gex.data[i,Her2.cols])
    luma.data <- as.numeric(gex.data[i,LumA.cols])
    lumb.data <- as.numeric(gex.data[i,LumB.cols])
    
    # for Her2 vs LumA
    res.A <- var_ttest(her2.data,luma.data)
    
    # save results
    Her2.LumA.pval[i] <- res.A$p.value
    Her2.LumA.diff[i] <- res.A$estimate[1]-res.A$estimate[2]
    
    # for Her2 vs LumB
    res.B <- var_ttest(her2.data,lumb.data)
    
    # save results
    Her2.LumB.pval[i] <- res.B$p.value
    Her2.LumB.diff[i] <- res.B$estimate[1]-res.B$estimate[2]
    close(pb)
}

# Process the output

# create the final output
DE.res <- gex.data %>% add_column(
    Her2.LumA.pval = Her2.LumA.pval,
    Her2.LumB.pval = Her2.LumB.pval,
    Her2.LumA.diff = Her2.LumA.diff, 
    Her2.LumB.diff = Her2.LumB.diff) %>% 
    dplyr::select(Her2.LumA.pval,Her2.LumB.pval,
                  Her2.LumA.diff,Her2.LumB.diff)

# adjust p value 
DE.res$Her2.LumA.padj <- p.adjust(DE.res$Her2.LumA.pval, "bonferroni")
DE.res$Her2.LumB.padj <- p.adjust(DE.res$Her2.LumB.pval, "bonferroni")

# up or down (refers to the situation in lumHer2 subtype; e.g. "up" indicates a gene whose expression is upregulated in lumHer2 compared to lumB or lumA)
DE.res <- DE.res %>% 
    mutate(Her2.LumA.de = case_when(
        Her2.LumA.diff <= 0 ~ "down", 
        Her2.LumA.diff >= 0 ~ "up")) %>% 
    mutate(Her2.LumB.de = case_when(
        Her2.LumB.diff <= 0 ~ "down", 
        Her2.LumB.diff >= 0 ~ "up"))

# save the file
save(DE.res,file = paste(data.path,"DE_results.RData",sep="") )
load(file = paste(data.path,"DE_results.RData",sep=""))