# Script: Generate UMAP

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
plot.file <- paste(output.path,cohort,"_HER2n_umap.pdf",sep = "")

#packages
source("scripts/2_transcriptomic/src/tscr_functions.R")
library(ggplot2)
library(tidyverse)
library(matrixStats)
library(VennDiagram)
library(readxl)
library(ggfortify)
library(janitor)
library(biomaRt)
library(amap)
library(umap)

#######################################################################
# 2. Cohort-specific data preprocessing including selection of  
# the clinical ER+Her2- subtyped samples
#######################################################################

# for SCANB
if (cohort=="SCANB") {
  
  # load annotation data and select subgroup data
  anno <- loadRData(file="./data/SCANB/1_clinical/processed/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data_ERpHER2n.RData") %>%
    dplyr::rename(sampleID = GEX.assay, PAM50 = NCN.PAM50) %>% 
    mutate(Point.size = ifelse(PAM50=="Her2",10,5)) # add point size for plotting
  
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
    dplyr::rename(sampleID=METABRIC_ID,NHG=Grade) %>% # rename to match SCANB variables
    mutate(Point.size = ifelse(PAM50=="Her2",10,5)) # add point size for plotting
  
  # load and select subgroup data
  gex.data <- metabric_gex_load("./data/METABRIC/2_transcriptomic/raw/data_mRNA_median_all_sample_Zscores.txt",ID.type = "Hugo_Symbol") %>% 
    dplyr::select(any_of(anno$sampleID)) %>% 
    mutate_all(function(x) as.numeric(x)) %>% 
    select_if(~ !any(is.na(.))) 
  
  # exclude samples from anno without associated gex data
  anno <- anno %>% 
    filter(sampleID %in% colnames(gex.data))
}

#######################################################################
# Load DE results 
#######################################################################

# load DE results
load(file = "./data/SCANB/2_transcriptomic/processed/DE_results.RData")
DE.res.scanb <- DE.res
load(file = "./data/METABRIC/2_transcriptomic/processed/DE_results.RData")
DE.res.metabric <- DE.res

#######################################################################
# UMAP with all core DEGs
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
deg.gex <- gex.data %>% 
  rownames_to_column("Gene") %>% 
  filter(Gene %in% DEGs.core) %>% 
  column_to_rownames(var="Gene") 

#---------------------------------------------------------------------#

# merge to one df to have pam50 annotation for each sample in the right order
deg.gex <- as.data.frame(t(deg.gex)) %>% 
  rownames_to_column(var="sampleID")
umap.gex <- merge(deg.gex,anno[c("PAM50","sampleID","Point.size")],by="sampleID") %>% 
  column_to_rownames(var="sampleID") 

umap.gex$PAM50 <- factor(umap.gex$PAM50,levels = c("LumA","LumB","Her2"))

# run umap
umap <- umap(umap.gex[1:(ncol(umap.gex)-2)],n_neighbors=50)
df <- data.frame(x = umap$layout[,1],
                 y = umap$layout[,2],
                 PAM50 = umap.gex["PAM50"],
                 Point.size = umap.gex["Point.size"]) %>% 
  arrange(PAM50)

# plot
plot <- ggplot(df, aes(x, y, colour = PAM50)) + #,order=
  geom_point(alpha=0.8, size=df$Point.size) + #
  theme_bw() +
  theme(axis.text.x = element_text(size = 40,margin = margin(t=10)),
        axis.title.x = element_text(size = 40),
        axis.text.y = element_text(size = 35,margin = margin(r=10)),
        axis.title.y = element_text(size = 40),
        plot.title = element_text(size=40),
        legend.position = "none",
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",linewidth=2),
        axis.ticks = element_line(colour = "black", linewidth = 2),
        axis.ticks.length=unit(0.5, "cm")) +
  scale_color_manual(values=setNames(c("#d334eb","#2176d5","#34c6eb"),
                                     c("Her2","LumA","LumB"))) +
  ggtitle(paste("UMAP based on core DEGs (n genes=",ncol(umap.gex)-2,"; ",cohort,")",sep=""))

# plot
plot.list <- append(plot.list, list(plot))

#######################################################################
# UMAP with all genes
#######################################################################

# merge to one df to have pam50 annotation for each sample in the right order
all.gex <- as.data.frame(t(gex.data)) %>% 
  rownames_to_column(var="sampleID")

umap.gex <- merge(all.gex,anno[c("PAM50","sampleID","Point.size")],by="sampleID") %>% 
  column_to_rownames(var="sampleID") 

umap.gex$PAM50 <- factor(umap.gex$PAM50,levels = c("LumA","LumB","Her2"))

# run umap
umap <- umap(umap.gex[1:(ncol(umap.gex)-2)],n_neighbors=50)

## export source dat
umap.sourcefile <- umap.gex[c("PAM50")]
umap.sourcefile$sampleID <- rownames(umap.sourcefile)
rownames(umap.sourcefile) <- NULL
umap.sourcefile <- umap.sourcefile[c("sampleID","PAM50")]
umap.sourcefile$UMAP1 <- umap$layout[,1]
umap.sourcefile$UMAP2 <- umap$layout[,2]
save(umap.sourcefile, file="./output/source_data/R_objects/Figure_3_umap.RData")
##

df <- data.frame(x = umap$layout[,1],
                 y = umap$layout[,2],
                 PAM50 = umap.gex["PAM50"],
                 Point.size = umap.gex["Point.size"]) %>% 
  arrange(PAM50)

# plot
plot <- ggplot(df, aes(x, y, colour = PAM50)) + #,order=
  geom_point(alpha=0.7, size=df$Point.size) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 40,margin = margin(t=10)),
        axis.title.x = element_text(size = 40),
        axis.text.y = element_text(size = 35,margin = margin(r=10)),
        axis.title.y = element_text(size = 40),
        plot.title = element_text(size=40),
        legend.position = "none",
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",linewidth=2),
        axis.ticks = element_line(colour = "black", linewidth = 2),
        axis.ticks.length=unit(0.5, "cm")) +
  scale_color_manual(values=setNames(c("#d334eb","#2176d5","#34c6eb"),
                                     c("Her2","LumA","LumB"))) +
  ggtitle(paste("UMAP based on all genes (n genes=",ncol(umap.gex)-2,"; ",cohort,")",sep=""))

# plot
plot.list <- append(plot.list, list(plot))

#######################################################################

# save plots
pdf(file = plot.file, 
    onefile = TRUE, width = 14.8, height = 10.5) 
for (i in 1:length(plot.list)) {
  print(plot.list[[i]])
}
dev.off()
