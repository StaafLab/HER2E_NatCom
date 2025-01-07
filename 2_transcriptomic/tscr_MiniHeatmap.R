# Script: Generate a mini heatmap for the 258 core DEGs
#TODO: 

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "SCANB" # SCANB or METABRIC

# set/create output directory for plots
output.path <- "output/plots/2_transcriptomic/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/2_transcriptomic/processed/",sep="")
dir.create(data.path)

# plot
plot.list <- list() # object to store plots; note: if the output is not in string format use capture.output()
plot.file <- paste(output.path,cohort,"_HER2n_miniHM.pdf",sep = "")

#packages
source("scripts/2_transcriptomic/src/tscr_functions.R")
library(ggplot2)
library(tidyverse)
library(matrixStats)
library(pheatmap)
library(ComplexHeatmap)
#library(Hmisc)
library(VennDiagram)
library(readxl)
library(ggfortify)
library(janitor)
library(biomaRt)
library(amap)

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
  anno <- loadRData("data/METABRIC/1_clinical/processed/Merged_annotations_ERpHER2n.RData") %>% 
    dplyr::rename(sampleID=METABRIC_ID,NHG=Grade) # rename to match SCANB variables
  
  # load and select subgroup data
  gex.data <- metabric_gex_load("./data/METABRIC/2_transcriptomic/raw/data_mRNA_median_all_sample_Zscores.txt",ID.type = "Hugo_Symbol") %>% 
    dplyr::select(any_of(anno$sampleID)) %>% 
    mutate_all(function(x) as.numeric(x))
  
  # exclude samples from anno without associated gex data
  anno <- anno %>% 
    filter(sampleID %in% colnames(gex.data))
}

#######################################################################
# Load DE results and convert ids
#######################################################################

# load DE results
load(file = "./data/SCANB/2_transcriptomic/processed/DE_results.RData")
DE.res.scanb <- DE.res
load(file = "./data/METABRIC/2_transcriptomic/processed/DE_results.RData")
DE.res.metabric <- DE.res

#######################################################################
# Heatmap parameters
#######################################################################

n.tree <- 5  # nbr of branches/clusters 
nbr.c <- 2 # color palette size = max value
breaksList = seq(-2, nbr.c, by = 0.1) # this defines the color palette range
my.link.method <- "ward.D" 
my.dist.method <- "euclidean" #"correlation" "euclidean"

# colors
my_colour = list(
  Cluster=c("1" = "#e41a1c","2" = "#377eb8","3" = "#4daf4a","4" = "#984ea3","5" = "#ff7f00"),
  HER2E = c("0" = "#f0f0f0" , "1" = "#000000"),
  LumA = c("0" = "#f0f0f0" , "1" = "#000000"),
  LumB = c("0" = "#f0f0f0" , "1" = "#000000"))
  #PAM50 = c("Her2" = "#d334eb" , "LumA" = "#2176d5" , "LumB" = "#34c6eb"))

#######################################################################
# Heatmap: Core DEGs
#######################################################################

# define DEG set
DEGs.scanb <- DE.res.scanb %>% 
  filter(Her2.LumA.padj <= 0.05) %>% 
  filter(Her2.LumB.padj <= 0.05) %>% 
  rownames_to_column("Gene") %>% pull(Gene)

DEGs.metabric <- DE.res.metabric %>% 
  filter(Her2.LumA.padj <= 0.05) %>% 
  filter(Her2.LumB.padj <= 0.05) %>% 
  rownames_to_column("Gene") %>% pull(Gene) 

# core gex
DEGs.core <- intersect(DEGs.scanb,DEGs.metabric)

# get gex of selected genes
core.gex <- gex.data %>% 
  rownames_to_column("Gene") %>% 
  filter(Gene %in% DEGs.core) %>% 
  column_to_rownames(var="Gene") 

# create hm annotaiton object
hm.anno <- anno[c("sampleID","PAM50")] %>% 
  column_to_rownames(var = "sampleID") %>% 
  # add extra columns for the seperate pam50 annotation in hm
  mutate(HER2E = ifelse(PAM50=="Her2",1,0)) %>% 
  mutate(LumA = ifelse(PAM50=="LumA",1,0)) %>%
  mutate(LumB = ifelse(PAM50=="LumB",1,0))

# cluster samples (prev. with dist())
c1 <- cutree(hclust(
  Dist(t(core.gex),method=my.dist.method),method=my.link.method),n.tree)

# add clusters to hm annotation object 
hm.anno <- merge(hm.anno,as.data.frame(c1),by=0) %>% 
  column_to_rownames(var = "Row.names") %>% 
  dplyr::rename(Cluster=c1)

# final hm annotation 
final.hm.anno <- data.frame(
  Cluster=factor(hm.anno$Cluster),
  PAM50=factor(hm.anno$PAM50),
  HER2E=factor(hm.anno$HER2E),
  LumA=factor(hm.anno$LumA),
  LumB=factor(hm.anno$LumB))
rownames(final.hm.anno) <- rownames(hm.anno)

# create the heatmap
plot1 <- pheatmap(core.gex, cluster_cols=T, cluster_rows=F, treeheight_row = 0,
                  clustering_distance_rows = my.dist.method,
                  clustering_distance_cols = my.dist.method, 
                  clustering_method = my.link.method,
                  cellwidth=0.3,
                  show_rownames=F, show_colnames=F, 
                  #main=paste("Core DEGs (ERpHER2n; ",cohort,")",sep = ""), 
                  cutree_cols=n.tree, legend = F,
                  #color=colorRampPalette(
                  # c("#998ec3", "#f7f7f7", "#f1a340"))(length(breaksList)),
                  annotation_col=final.hm.anno[rev(c("Cluster","HER2E","LumA","LumB"))], 
                  annotation_colors=my_colour,breaks=breaksList)

#print(plot1)
plot.list <- append(plot.list, list(plot1))

#######################################################################
#######################################################################

# save plots
pdf(file = plot.file, onefile = TRUE)

for (i in 1:length(plot.list)) {
  grid::grid.newpage()
  grid::grid.draw(plot.list[[i]]$gtable)
  # try to draw a white rectable over the irrelevant info
  #grid::grid.draw(rectGrob(x = 0, y = 0, width = 1, height = 0.9,gp=gpar(fill="black")))
  #print(plot.list[[i]])
}

dev.off()

