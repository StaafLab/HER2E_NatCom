# Script: Pathway enrichment analysis of PAM50 gene clusters

# TODO:

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "SCANB" # SCANB 

# set/create output directory for plots
output.path <- "output/plots/1_clinical/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/1_clinical/processed/",sep="")
dir.create(data.path)

#packages
#source("scripts/1_clinical/src/clin_functions.R")
library(ggplot2)
library(tidyverse)
library(reshape)
library(openxlsx)
library(org.Hs.eg.db)
library(enrichR)

#######################################################################
# Core clusters
#######################################################################

# load data
clust.anno <- as.data.frame(read_table("data/SCANB/1_clinical/raw/Core_Clusters_6_Info_ann.txt"))

# convert hgnc to entrez ids 
# hgnc.ids <- clust.anno %>% pull(HGNC) 
# hs <- org.Hs.eg.db
# res <- select(hs, 
#               keys = hgnc.ids,
#               columns = c("ENTREZID", "SYMBOL"),
#               keytype = "SYMBOL") %>% dplyr::rename(HGNC=SYMBOL)

# merge
#clust.anno <- as.data.frame(merge(clust.anno,res,by="HGNC"))

# pwe analysis
setEnrichrSite("Enrichr") 
dbs <- listEnrichrDbs() # see avaiable dbs and select 
dbs <- c("KEGG_2021_Human") #"GO_Molecular_Function_2021", "GO_Biological_Process_2021

# for each cluster
pdf(file = paste(output.path,cohort,"_PAM50clusters_pwe.pdf", sep=""), onefile = TRUE )
for(i in 1:length(unique(clust.anno$ClusterNumber))) { 
    data <- clust.anno %>% filter(ClusterNumber == i)
    clust <- paste("core.cluster.", i, sep = "")
    #assign(clust, enrichr(data$HGNC, dbs)) # will seperate for plotting
    res <- enrichr(data$HGNC, dbs)
    if (nrow(res[[1]])>0) {
        plot <- plotEnrich(res[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value", title = paste(clust))
        print(plot)
    } else {}
    assign(clust, res)
}

#######################################################################
# Spiral clusters
#######################################################################

# load data
clust.anno <- as.data.frame(read_table("data/SCANB/1_clinical/raw/Spiral_SRIQ_PAM50_Gex_Clusters_6_Info_ann.txt"))

# convert hgnc to entrez ids 
# hgnc.ids <- clust.anno %>% pull(HGNC) 
# hs <- org.Hs.eg.db
# res <- select(hs, 
#               keys = hgnc.ids,
#               columns = c("ENTREZID", "SYMBOL"),
#               keytype = "SYMBOL") %>% dplyr::rename(HGNC=SYMBOL)

# merge
#clust.anno <- as.data.frame(merge(clust.anno,res,by="HGNC"))

# pwe analysis
setEnrichrSite("Enrichr") 
dbs <- listEnrichrDbs() # see avaiable dbs and select 
dbs <- c("KEGG_2021_Human") #"GO_Molecular_Function_2021", "GO_Biological_Process_2021

# for each cluster
for(i in 1:length(unique(clust.anno$ClusterNumber))) { 
    data <- clust.anno %>% filter(ClusterNumber == i)
    clust <- paste("spiral.cluster.", i, sep = "")
    #assign(clust, enrichr(data$HGNC, dbs)) # will seperate for plotting
    res <- enrichr(data$HGNC, dbs)
    if (nrow(res[[1]])>0) {
        plot <- plotEnrich(res[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value", title = paste(clust))
        print(plot)
    } else {}
    assign(clust, res)
}

dev.off()
