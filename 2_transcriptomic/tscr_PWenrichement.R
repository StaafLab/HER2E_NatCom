# Script: Pathway enrichment analysis for DEGs

# TODO

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "COMBINED" 

# set/create output directory for plots
output.path <- "output/plots/2_transcriptomic/"
dir.create(output.path)

# output filenames
plot.list <- list() # object to store plots; note: if the output is not in string format use capture.output()
plot.file <- paste(output.path,cohort,"_HER2n_pwenrichment_plots_wikipathways.pdf", sep="")
xlsx.list <- list() # object to store dataframes -> excel tables 
xlsx.file <- paste("output/supplementary_data/HER2n_pwenrichment_wikipathways.xlsx", sep="")

#packages
source("./scripts/2_transcriptomic/src/tscr_functions.R")
library(ggplot2)
library(tidyverse)
library(reshape)
library(openxlsx)
library(enrichR)
library(janitor)

#######################################################################

# load data
DE.res.scanb <- loadRData("./data/SCANB/2_transcriptomic/processed/DE_results.RData")
DE.res.metabric <- loadRData("./data/METABRIC/2_transcriptomic/processed/DE_results.RData")

# set database parameters
setEnrichrSite("Enrichr") 
dbs <- listEnrichrDbs() # see avaiable dbs and select 
dbs <- c("WikiPathway_2023_Human") #"GO_Biological_Process_2021, "WikiPathway_2023_Human"
plot.file <- paste(output.path,cohort,"_HER2n_pwenrichment_plots_",dbs,".pdf", sep="")
xlsx.list <- list() # object to store dataframes -> excel tables 
xlsx.file <- paste("output/supplementary_data/HER2n_pwenrichment_",dbs,".xlsx", sep="")

#######################################################################
# Core DEGs
#######################################################################

#filter based on expression diff?
#filter(abs(Her2.LumA.diff) >= 1) %>% 

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

#---------------------------------------------------------------------#

#assign("res", enrichr(top.DEGs.core, dbs)) 

res.core <- as.data.frame(enrichr(DEGs.core, dbs)[[1]])
xlsx.list <- append(xlsx.list, list("CoreDEGs" = res.core))

#######################################################################
# HER2E vs. LumAs
#######################################################################

# define DEG set
LumA.DEGs.scanb <- DE.res.scanb %>% 
  filter(Her2.LumA.padj <= 0.05) %>% 
  rownames_to_column("Gene") %>% 
  pull(Gene)

LumA.DEGs.metabric <- DE.res.metabric %>% 
  filter(Her2.LumA.padj <= 0.05) %>% 
  rownames_to_column("Gene") %>% 
  pull(Gene) 

# gex
LumA.DEGs.core <- intersect(LumA.DEGs.scanb,LumA.DEGs.metabric)

#---------------------------------------------------------------------#

res.luma <- as.data.frame(enrichr(LumA.DEGs.core, dbs)[[1]])
xlsx.list <- append(xlsx.list, list("LumA_coreDEGs" = res.luma))

#######################################################################
# HER2E vs. LumB
#######################################################################

# define DEG set
LumB.DEGs.scanb <- DE.res.scanb %>% 
  filter(Her2.LumB.padj <= 0.05) %>% 
  rownames_to_column("Gene") %>% 
  pull(Gene)

LumB.DEGs.metabric <- DE.res.metabric %>% 
  filter(Her2.LumB.padj <= 0.05) %>% 
  rownames_to_column("Gene") %>% 
  pull(Gene) 

# core top gex
LumB.DEGs.core <- intersect(LumB.DEGs.scanb,LumB.DEGs.metabric)

#---------------------------------------------------------------------#

res.lumb <- as.data.frame(enrichr(LumB.DEGs.core, dbs)[[1]])
xlsx.list <- append(xlsx.list, list("LumB_coreDEGs" = res.lumb))

#######################################################################
# plotting
#######################################################################

# filter and get ready for combining into single plot file
#core
res.core <- res.core[res.core$Adjusted.P.value <= 0.05,] %>% # filter signif
  mutate(Comp="Core") %>%  # add column for group
  mutate(Gene_count = as.numeric(sapply(strsplit(Overlap, "/"), "[[", 1))) %>% # get gene count from overlap column
  arrange(desc(Gene_count),Adjusted.P.value) %>% # order
  mutate(Adjusted.P.value= round(Adjusted.P.value,3)) %>% # round pval
  slice(1:3) # select top n rows

#luma
res.luma <- res.luma[res.luma$Adjusted.P.value <= 0.05,] %>% 
  mutate(Comp="LumA") %>% 
  mutate(Gene_count = as.numeric(sapply(strsplit(Overlap, "/"), "[[", 1))) %>% 
  arrange(desc(Gene_count),Adjusted.P.value) %>% 
  mutate(Adjusted.P.value= round(Adjusted.P.value,3)) %>% 
  slice(1:3)

#lumb
res.lumb <- res.lumb[res.lumb$Adjusted.P.value <= 0.05,] %>% 
  mutate(Comp="LumB") %>% 
  mutate(Gene_count = as.numeric(sapply(strsplit(Overlap, "/"), "[[", 1))) %>% 
  arrange(desc(Gene_count),Adjusted.P.value) %>% 
  mutate(Adjusted.P.value= round(Adjusted.P.value,3)) %>% 
  slice(1:3)

res.all <- rbind(res.core,res.luma,res.lumb) %>% 
  mutate(Term = gsub("\\s*\\([^\\)]+\\)","",as.character(Term)))
res.all$Comp <- as.factor(res.all$Comp)
#View(res.all)

# labeller function for plot
deg.sets <- list("Core" = paste("Core DEGs (n=",length(DEGs.core),")",sep=""),
                 "LumA" = paste("LumA DEGs (n=",length(LumA.DEGs.core),")",sep=""),
                 "LumB" = paste("LumB DEGs (n=",length(LumB.DEGs.core),")",sep=""))

plot.list <- append(plot.list, list(
  pwplot(res.all,title="Enriched pathways in DEG sets")))

#######################################################################
# top core DEGs
#######################################################################
dbs <- c("WikiPathway_2023_Human") 
tc.degs <- c("ESR1","TBC1D9","CCDC170","RERG","IGF1R","FGFR4","BCL2","SOX11","THSD4")

head(as.data.frame(enrichr(tc.degs, dbs)[[1]]))

#######################################################################
# save plots and excel file
#######################################################################

# save excel file
write.xlsx(xlsx.list, file = xlsx.file)

# save plots
pdf(file = plot.file, onefile = TRUE, width = 25, height = 14) #

for (i in 1:length(plot.list)) {
  print(plot.list[[i]])
}

dev.off()
