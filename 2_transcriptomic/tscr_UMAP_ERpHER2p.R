# Script: Generate UMAP in HER2-positive

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

# plot
plot.list <- list() # object to store plots; note: if the output is not in string format use capture.output()
plot.file <- paste(output.path,cohort,"_HER2p_umap.pdf",sep = "")

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

######################################################################
# Cohort-specific data preprocessing 
#######################################################################

# for SCANB
if (cohort=="SCANB") {
  
  # load annotation data and select subgroup data
  anno <- loadRData("./data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData") %>% 
    filter(Follow.up.cohort==TRUE) %>% 
    filter(fuV8==TRUE) %>% 
    filter(ER=="Positive") %>% 
    dplyr::rename(sampleID = GEX.assay, PAM50 = NCN.PAM50) %>% 
    mutate(Group = case_when(
      HER2 == "Negative" & PAM50 == "Her2" ~ "HER2n_HER2E",
      HER2 == "Positive" & PAM50 == "Her2" ~ "HER2p_HER2E",
      HER2 == "Positive" & PAM50 != "Her2" ~ "HER2p_nonHER2E")) %>% 
    filter(Group %in% 
             c("HER2n_HER2E","HER2p_HER2E","HER2p_nonHER2E")) %>% 
    mutate(Point.size = ifelse(Group=="HER2n_HER2E",10,5)) # add point size for plotting
  
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
    mutate(HER2 = case_when(
      HER2_SNP6_state == "GAIN" ~ "Positive",
      HER2_SNP6_state == "UNDEF" ~ "Undefined",
      HER2_SNP6_state %in% c("LOSS","NEUT") ~ "Negative")) %>%
    filter(ER == "Positive") %>% 
    mutate(Group = case_when(
      HER2 == "Negative" & PAM50 == "Her2" ~ "HER2n_HER2E",
      HER2 == "Positive" & PAM50 == "Her2" ~ "HER2p_HER2E",
      HER2 == "Positive" & PAM50 != "Her2" ~ "HER2p_nonHER2E")) %>% 
    filter(Group %in% 
             c("HER2n_HER2E","HER2p_HER2E","HER2p_nonHER2E")) %>% 
    dplyr::rename(sampleID=METABRIC_ID)  %>%  #rename to match SCANB variables
    mutate(Point.size = ifelse(Group=="HER2n_HER2E",10,5)) # add point size for plotting#
  
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
# UMAP with all genes
#######################################################################

# merge to one df to have pam50 annotation for each sample in the right order
all.gex <- as.data.frame(t(gex.data)) %>% 
  rownames_to_column(var="sampleID")

umap.gex <- merge(all.gex,anno[c("Group","sampleID","Point.size")],by="sampleID") %>% 
  column_to_rownames(var="sampleID") 

umap.gex$Group <- factor(umap.gex$Group,levels = c("HER2p_nonHER2E", 
                                                   "HER2p_HER2E","HER2n_HER2E"))

# run umap
umap <- umap(umap.gex[1:(ncol(umap.gex)-2)],n_neighbors=50) 

## export source dat
umap.sourcefile <- umap.gex[c("Group")]
umap.sourcefile$sampleID <- rownames(umap.sourcefile)
rownames(umap.sourcefile) <- NULL
umap.sourcefile <- umap.sourcefile[c("sampleID","Group")]
umap.sourcefile$UMAP1 <- umap$layout[,1]
umap.sourcefile$UMAP2 <- umap$layout[,2]
save(umap.sourcefile, file="./output/source_data/R_objects/Figure_5_umapHER2p.RData")
##



df <- data.frame(x = umap$layout[,1],
                 y = umap$layout[,2],
                 Group = umap.gex["Group"],
                 Point.size = umap.gex["Point.size"]) %>% 
  arrange(Group)

# plot
plot <- ggplot(df, aes(x, y, colour = Group)) + #,order=
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
  scale_color_manual(values=setNames(c("#d334eb","#d8b365","#5ab4ac"),
                                     c("HER2n_HER2E","HER2p_nonHER2E", 
                                       "HER2p_HER2E"))) +
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
