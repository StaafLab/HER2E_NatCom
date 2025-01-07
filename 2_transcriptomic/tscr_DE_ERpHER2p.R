# Script: differential expression analysis in HER2p ; cohorts=SCANB, METABRIC
#TODO: 

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
    anno <- as.data.frame(
        read_excel("data/SCANB/1_clinical/raw/NPJ_release.xlsx")) %>%
        filter(Follow.up.cohort==TRUE) %>% 
        filter(ER=="Positive") %>% 
        dplyr::rename(sampleID = GEX.assay, PAM50 = NCN.PAM50) %>% 
        mutate(Group = case_when(
            HER2 == "Negative" & PAM50 == "Her2" ~ "HER2n_HER2E",
            HER2 == "Positive" & PAM50 == "Her2" ~ "HER2p_HER2E",
            HER2 == "Positive" & PAM50 != "Her2" ~ "HER2p_nonHER2E")) %>% 
        filter(Group %in% 
                   c("HER2n_HER2E","HER2p_HER2E","HER2p_nonHER2E"))
    
    # load gex data
    gex.data <- scanb_gex_load(gex.path = "data/SCANB/2_transcriptomic/raw/genematrix_noNeg.Rdata", geneanno.path = "data/SCANB/1_clinical/raw/Gene.ID.ann.Rdata", ID.type = "Gene.Name") %>% 
        dplyr::select(any_of(anno$sampleID)) %>% # select subgroup gex 
        select_if(~ !any(is.na(.))) # otherwise error when scaling
    
    # log transform FPKM data
    gex.data <- as.data.frame(log2(gex.data + 1))
    # z-transform
    gex.data <- as.data.frame(t(apply(gex.data, 1, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y))))) # for some rows there may be 0 variance so i have to handle these cases
    
#-----------------------------------------------------------------------#
    
} else if (cohort=="METABRIC") {
    
    # load annotation data
    load("data/METABRIC/1_clinical/raw/Merged_annotations.RData")
    
    # extract relevant variables
    anno <- anno %>% 
        mutate(ER = if_else(
            ER_IHC_status=="pos","Positive","Negative")) %>% 
        mutate(HER2 = if_else(
            HER2_SNP6_state=="GAIN","Positive","Negative")) %>% 
        filter(ER == "Positive") %>% 
        mutate(Group = case_when(
            HER2 == "Negative" & PAM50 == "Her2" ~ "HER2n_HER2E",
            HER2 == "Positive" & PAM50 == "Her2" ~ "HER2p_HER2E",
            HER2 == "Positive" & PAM50 != "Her2" ~ "HER2p_nonHER2E")) %>% 
        filter(Group %in% 
                   c("HER2n_HER2E","HER2p_HER2E","HER2p_nonHER2E")) %>% 
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
g1.cols <- anno %>% filter(Group=="HER2n_HER2E") %>% pull(sampleID) #HER2n_HER2E
g2.cols <- anno %>% filter(Group=="HER2p_nonHER2E") %>% pull(sampleID) # luma HER2p_nonHER2E
g3.cols <- anno %>% filter(Group=="HER2p_HER2E") %>% pull(sampleID) #lumb g3 HER2p_HER2E

# initialize vector of stored p-values and expression differences (based on mean comparison)
g1.g2.pval <- rep(NA,nrow(gex.data))
g1.g3.pval <- rep(NA,nrow(gex.data))
g1.g2.diff <- rep(NA,nrow(gex.data))
g1.g3.diff <- rep(NA,nrow(gex.data))

# loop through the genes
pb = txtProgressBar(min = 0, max = nrow(gex.data), initial = 0, style = 3)
for (i in 1:nrow(gex.data)) { #nrow(gex.data)
    setTxtProgressBar(pb,i)
    
    # set vars
    g1.data <- as.numeric(gex.data[i,g1.cols])
    g2.data <- as.numeric(gex.data[i,g2.cols])
    g3.data <- as.numeric(gex.data[i,g3.cols])
    
    # for g1 vs g2
    res.A <- var_ttest(g1.data,g2.data)
    
    # save results
    g1.g2.pval[i] <- res.A$p.value
    g1.g2.diff[i] <- res.A$estimate[1]-res.A$estimate[2]
    
    # for g1 vs g3
    res.B <- var_ttest(g1.data,g3.data)
    
    # save results
    g1.g3.pval[i] <- res.B$p.value
    g1.g3.diff[i] <- res.B$estimate[1]-res.B$estimate[2]
    close(pb)
}

# Process the output

# create the final output
DE.res <- gex.data %>% add_column(
    g1.g2.pval = g1.g2.pval,
    g1.g3.pval = g1.g3.pval,
    g1.g2.diff = g1.g2.diff, 
    g1.g3.diff = g1.g3.diff) %>% 
    dplyr::select(g1.g2.pval,g1.g3.pval,
                  g1.g2.diff,g1.g3.diff)

# adjust p value 
DE.res$g1.g2.padj <- p.adjust(DE.res$g1.g2.pval, "bonferroni")
DE.res$g1.g3.padj <- p.adjust(DE.res$g1.g3.pval, "bonferroni")

# up or down (refers to the situation in lumHer2 subtype; e.g. "up" indicates a gene whose expression is upregulated in lumHer2 compared to lumB or lumA)
DE.res <- DE.res %>% 
    mutate(g1.g2.de = case_when(
        g1.g2.diff <= 0 ~ "down", 
        g1.g2.diff >= 0 ~ "up")) %>% 
    mutate(g1.g3.de = case_when(
        g1.g3.diff <= 0 ~ "down", 
        g1.g3.diff >= 0 ~ "up"))

# rename columns
colnames(DE.res) <- gsub("g1", "HER2n_HER2E", colnames(DE.res))
colnames(DE.res) <- gsub("g2", "HER2p_nonHER2E", colnames(DE.res))
colnames(DE.res) <- gsub("g3", "HER2p_HER2E", colnames(DE.res))

# save the file
save(DE.res,file = paste(data.path,"HER2p_DE_results.RData",sep="") )
#load(file = paste(data.path,"DE_results.RData",sep=""))
