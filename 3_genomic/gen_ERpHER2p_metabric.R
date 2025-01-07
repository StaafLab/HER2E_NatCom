# Script: Mutational analyses in Metabric (not WGS-based)

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "METABRIC"

# set/create output directory for plots
output.path <- "output/plots/3_genomic/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/3_genomic/processed/",sep="")
dir.create(data.path)

# plot
plot.list <- list() # object to store plots; note: if the output is not in string format use capture.output()
plot.file <- paste("output/plots/3_genomic/METABRIC_HER2p_mutbarplots.pdf",sep = "")
txt.out <- c() # object to store text output
txt.file <- paste(output.path,cohort,"_HER2p_mutbarplots.txt", sep="")

#packages
source("scripts/3_genomic/src/gen_functions.R")
source("scripts/4_CN/src/cn_functions.R")
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
library(ggstatsplot)
library(GenVisR)

#######################################################################
# functions
#######################################################################

# function to convert the freq table to long format
countsToCases <- function(x, countcol = "Freq") {
    # Get the row indices to pull from x
    idx <- rep.int(seq_len(nrow(x)), x[[countcol]])
    # Drop count column
    x[[countcol]] <- NULL
    # Get the rows from x
    x[idx, ]
}

# 1 mut per sample
count.sample <- function(data,gene,group.n) {
  data <- data %>% 
    distinct(Sample,Gene, .keep_all = TRUE)
  data[data$Gene==gene,] %>% 
    dplyr::count(Group) %>% 
    mutate(Freq=case_when(Group=="HER2n_HER2E" ~ (n/group.n["HER2n_HER2E"])*100,
                          Group=="HER2p_nonHER2E" ~ (n/group.n["HER2p_nonHER2E"])*100,
                          Group=="HER2p_HER2E" ~ (n/group.n["HER2p_HER2E"])*100))
}

#######################################################################
# add anno column
#######################################################################

# # add her2 amp status to anno
# # extract relevant variables
# anno <- loadRData("data/METABRIC/1_clinical/raw/Merged_annotations.RData")
# cn.data <- read.delim("./data/METABRIC/4_CN/raw/data_CNA.txt",header = TRUE, sep = "\t", dec = ".") %>% 
#   dplyr::select(-c(Entrez_Gene_Id)) 
# names(cn.data) <- gsub(x = names(cn.data), 
#                        pattern = "\\.", 
#                        replacement = "-") 
# cn.data <- cn.data %>% 
#   filter(Hugo_Symbol=="ERBB2") %>% 
#   dplyr::select(any_of(c("Hugo_Symbol",anno$METABRIC_ID))) %>% 
#   pivot_longer(cols= any_of(anno$METABRIC_ID), 
#                names_to = "METABRIC_ID", 
#                values_to="HER2_amp") %>% 
#   dplyr::select(-c("Hugo_Symbol")) %>% 
#   mutate(HER2_amp = ifelse(HER2_amp==2,"yes","no"))
# # join on sample name sample HER2_amp
# anno <- as.data.frame(merge(anno,cn.data,by="METABRIC_ID",all=TRUE))
# save(anno, file= "data/METABRIC/1_clinical/raw/Merged_annotations.RData")

#######################################################################
# 2. load data
#######################################################################

# genes 
driver.genes <- loadRData("data/SCANB/3_genomic/processed/driver_mutations_all.RData") %>% pull(gene) %>% unique()

# load data
mut.data <- as.data.frame(read.delim('data/METABRIC/3_genomic/raw/data_mutations_extended.txt', header = FALSE, sep = "\t", dec = "."))

mut.data <- mut.data[-1,] %>% 
  row_to_names(row_number = 1) %>% 
  dplyr::rename(gene=Hugo_Symbol,variant_class=Variant_Classification,
                sample=Tumor_Sample_Barcode) %>% 
  dplyr::select(sample,gene,variant_class)

# extract relevant variables
anno <- loadRData("data/METABRIC/1_clinical/raw/Merged_annotations.RData") %>% 
  mutate(ER = if_else(
    ER_IHC_status=="pos","Positive","Negative")) %>% 
  filter(ER == "Positive") %>% 
  mutate(Group = case_when(
    ClinGroup == "ERpHER2n" & PAM50 == "Her2" ~ "HER2n_HER2E",
    HER2_amp == "yes" & PAM50 == "Her2" ~ "HER2p_HER2E",
    HER2_amp == "yes" & PAM50 != "Her2" ~ "HER2p_nonHER2E")) %>% 
  filter(Group %in% 
           c("HER2n_HER2E","HER2p_HER2E","HER2p_nonHER2E")) %>% 
  dplyr::rename(sample=METABRIC_ID) %>% # rename to match SCANB variables
  dplyr::select(sample,Group)
    
mut.data <- mut.data %>% filter(sample %in% anno$sample)

# add amp alterations to this 
cn.data <- read.delim("./data/METABRIC/4_CN/raw/data_CNA.txt",header = TRUE, sep = "\t", dec = ".") %>% 
  dplyr::select(-c(Entrez_Gene_Id)) 
names(cn.data) <- gsub(x = names(cn.data), 
                       pattern = "\\.", 
                       replacement = "-") 
cn.data <- cn.data %>% 
  dplyr::select(any_of(c("Hugo_Symbol",anno$sample))) %>% 
  pivot_longer(cols= any_of(anno$sample), names_to = "sample", values_to="variant_class") %>% 
  dplyr::rename("gene"="Hugo_Symbol") %>% 
  filter(variant_class == 2) %>% 
  mutate(variant_class = "amplified")

mut.data <- rbind(mut.data,cn.data) %>% filter(gene %in% driver.genes)

#######################################################################
# 3. Waterfall plots
#######################################################################

# colors <- c("#8dd3c7",
#             "#ffffb3",
#             "#bebada",
#             "#fb8072",
#             "#80b1d3",
#             "#fdb462",
#             "#b3de69",
#             "#fccde5")

# Create a vector to save mutation priority order for plotting
mutation.priority <- as.character(unique(mut.data$variant_class))

mutation.priority <- c("amplified","Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","Splice_Site","In_Frame_Del","In_Frame_Ins","Splice_Region","Intron","3'UTR","Nonstop_Mutation","Translation_Start_Site","3'Flank","Silent")

# her2 wf plot
Her2p.samples <- anno %>% filter(Group=="HER2p_HER2E") %>% pull(sample)
mut.data.her2p <- mut.data %>% filter(sample %in% Her2p.samples)

## source dat exprot
save(mut.data.her2p,file="./output/source_data/R_objects/Figure_5_wfHER2p.RData")
##

pdf(file = "output/plots/3_genomic/METABRIC_HER2pHER2E_WF.pdf", onefile = TRUE,height = 10, width = 20)#, height = 10, width = 15)


wf.plot <- waterfall(mut.data.her2p, 
                     
                     fileType = "Custom", 
                     variant_class_order = mutation.priority,
                     mainGrid = TRUE,
                     mainPalette = c("#8dd3c7",
                                     "#ffffb3",
                                     "#bebada",
                                     "#fb8072",
                                     "#80b1d3",
                                     "#fdb462",
                                     "#b3de69",
                                     "#fccde5",
                                     "#a432a8"),
                     main_geneLabSize = 15,
                     mainRecurCutoff = 0,
                     maxGenes = 5,
                     mainDropMut = TRUE, # drop unused mutation types from legend
                     #rmvSilent = TRUE,
                     plotSamples = anno[anno$Group=="HER2p_HER2E",]$sample,
                     out= "grob",
                     mutBurdenLayer = layer,
                     plotMutBurden = FALSE)

grid::grid.newpage()
grid::grid.draw(wf.plot)

dev.off()

#######################################################################
# plot genes of interest: mutation frequencies
#######################################################################

# annotate samples 
mut.data$Group <- anno$Group[match(mut.data$sample,anno$sample)]

# adapt colnames 
mut.data <- mut.data %>% 
  dplyr::rename(Sample=sample,Gene=gene)

## source file export
# mut.data.sf <- mut.data
# mut.data.sf <- mut.data.sf[which(mut.data.sf$Gene %in% c("TP53","ERBB2")),]
# mut.data.sf <- mut.data.sf[!(mut.data.sf$Gene=="ERBB2" & mut.data.sf$variant_class == "amplified"),]
# mut.data.sf <- mut.data.sf %>% 
#   distinct(Sample,Gene, .keep_all = TRUE)
# mut.data.sf$variant_class <- NULL
# 
# # Create a new data frame with unique samples
# mut.data.sf.binary <- unique(mut.data.sf[, c("Sample", "Group")])
# 
# # Create binary columns for ERBB2 and TP53 mutations
# mut.data.sf.binary$ERBB2mut <- as.integer(mut.data.sf.binary$Sample %in% mut.data.sf$Sample[mut.data.sf$Gene == "ERBB2"])
# mut.data.sf.binary$TP53mut  <- as.integer(mut.data.sf.binary$Sample %in% mut.data.sf$Sample[mut.data.sf$Gene == "TP53"])
# 
# anno.sf <- anno
# colnames(anno.sf)[colnames(anno.sf) == "sample"] <- "Sample"
# # Merge mut.data.sf.binary and anno
# mut.data.sf.binary.all <- merge(anno.sf, mut.data.sf.binary, by = c("Sample", "Group"), all.x = TRUE)

# Replace NA values in ERBB2mut and TP53mut with 0 for samples not in mut.data.sf.binary
#mut.data.sf.binary.all$ERBB2mut[is.na(mut.data.sf.binary.all$ERBB2mut)] <- 0
#mut.data.sf.binary.all$TP53mut[is.na(mut.data.sf.binary.all$TP53mut)] <- 0
#save(mut.data.sf.binary.all, file="./output/source_data/R_objects/Figure_5_mutfreqsHER2p.RData")
#total_samples <- table(mut.data.sf.binary.all$Group)
#tapply(mut.data.sf.binary.all$ERBB2mut, mut.data.sf.binary.all$Group, mean)*100
#tapply(mut.data.sf.binary.all$TP53mut, mut.data.sf.binary.all$Group, mean)*100

##

# sample numbers
group.n <- table(
  mut.data[!duplicated(mut.data[,c("Sample")]),]$Group) 

# selected genes
gene.vec <-c("ERBB2","TP53","PIK3CA","MYC")

for (g in gene.vec) {
  
  if (g=="ERBB2") { # dont ocunt amp alterations duh
    plot.data <- count.sample(mut.data[which(mut.data$variant_class != "amplified"),],gene=g,group.n)
  } else {
    plot.data <- count.sample(mut.data,gene=g,group.n)
  }
  
  # Plots freq samples mutated in Group groups
  p2 <- ggplot(plot.data, aes(fill=as.factor(Group),x=Group,y=Freq)) +
    geom_bar(position="stack", stat="identity") +
    ggtitle(g) +
    theme_bw() +
    theme(aspect.ratio=1/1,
          legend.position = "none",
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black",linewidth = 2),
          axis.ticks = element_line(colour = "black", linewidth = 1),
          axis.ticks.length=unit(0.2, "cm")) +
    scale_fill_manual(values=setNames(c("#d334eb","#d8b365","#5ab4ac"),
                                      c("HER2n_HER2E","HER2p_nonHER2E", 
                                        "HER2p_HER2E"))) +
    scale_y_continuous(breaks = scales::breaks_pretty(10),
                       limits = c(0,80)) +
    ylab("Mutation frequency (%)") +
    xlab("Group subtype")
  
  plot.list <- append(plot.list,list(p2)) 
  txt.out <- append(txt.out, c(g))
  txt.out <- append(txt.out, c(capture.output(plot.data)))
}

#######################################################################
# stat testing mut. freqs
#######################################################################

# need df in binary format with samples as columns and genes as rows
sample.anno <- mut.data[c("Sample","Group")] %>% distinct(Sample, .keep_all = TRUE)
all.dmut.binary <- mut.data %>% 
  dplyr::select(-c(variant_class, Group)) %>% 
  distinct(Sample,Gene, .keep_all = TRUE) %>% 
  mutate(Mutation = 1) %>%  
  pivot_wider(names_from = Sample, values_from = Mutation) %>% 
  replace(is.na(.), 0) %>% 
  filter(Gene %in% gene.vec) %>% 
  filter(Gene != "ERBB2") # remove erbb2 to replce in next step

# replace erbb2 row without the amp alterations
erbb2.row <- mut.data[which(mut.data$Gene == "ERBB2" & mut.data$variant_class != "amplified"),c("Sample","Gene")] %>% 
  distinct(Sample,Gene, .keep_all = TRUE) %>% 
  mutate(Mutation = 1) %>%  
  pivot_wider(names_from = Sample, values_from = Mutation)

all.dmut.binary <- bind_rows(all.dmut.binary, erbb2.row) %>% 
  replace(is.na(.), 0)

res.df <- data.frame()

pb = txtProgressBar(min = 0, max = nrow(all.dmut.binary), initial = 0, style = 3)
for (i in 1:nrow(all.dmut.binary)) { #nrow(gex.data)
  setTxtProgressBar(pb,i)
  
  # gene to test
  gene <- all.dmut.binary$Gene[i]
  
  # mutation counts
  HER2n_HER2E.n.mut <- sum(all.dmut.binary[all.dmut.binary$Gene==gene,
                  sample.anno[sample.anno$Group=="HER2n_HER2E","Sample"]]==1)
  
  HER2p_nonHER2E.n.mut <- sum(all.dmut.binary[all.dmut.binary$Gene==gene,
                  sample.anno[sample.anno$Group=="HER2p_nonHER2E","Sample"]]==1)
  
  HER2p_HER2E.n.mut <- sum(all.dmut.binary[all.dmut.binary$Gene==gene,
                  sample.anno[sample.anno$Group=="HER2p_HER2E","Sample"]]==1)
  
  # make tbl
  freq.tbl <- data.frame(
    HER2n_HER2E = 
      c(HER2n_HER2E.n.mut,length(sample.anno[sample.anno$Group=="HER2n_HER2E","Sample"])-HER2n_HER2E.n.mut),
    HER2p_nonHER2E = 
      c(HER2p_nonHER2E.n.mut,length(sample.anno[sample.anno$Group=="HER2p_nonHER2E","Sample"])-HER2p_nonHER2E.n.mut),
    HER2p_HER2E = 
      c(HER2p_HER2E.n.mut,length(sample.anno[sample.anno$Group=="HER2p_HER2E","Sample"])-HER2p_HER2E.n.mut),
    row.names = c("mutation", "no_mutation"))
  
  # HER2n_HER2E vs HER2p_nonHER2E
  HER2p_nonHER2E.freq.tbl <- freq.tbl[,c("HER2n_HER2E","HER2p_nonHER2E")]
  HER2p_nonHER2E.res <- fisher.test(HER2p_nonHER2E.freq.tbl)
  
  # for HER2n_HER2E vs HER2p_HER2E
  HER2p_HER2E.freq.tbl <- freq.tbl[,c("HER2n_HER2E","HER2p_HER2E")]
  HER2p_HER2E.res <- fisher.test(HER2p_HER2E.freq.tbl)
  
  # save results 
  res.df <- rbind(res.df,c(
    gene,HER2p_nonHER2E.res$p.value,HER2p_HER2E.res$p.value))
  txt.out <- append(txt.out,c(gene,capture.output(HER2p_nonHER2E.res)))
  txt.out <- append(txt.out,c(gene,capture.output(HER2p_HER2E.res)))
  close(pb)
}

# name output columns
names(res.df) <- c("gene","HER2p_nonHER2E.pval","HER2p_HER2E.pval")

#######################################################################
#######################################################################

# save plots
pdf(file = plot.file, onefile = TRUE)#, height = 10, width = 15)

for (i in 1:length(plot.list)) {
  grid::grid.newpage()
  grid::grid.draw(plot.list[[i]])
  
  #print(plot.list[[i]])
}

dev.off()
# save text output
writeLines(txt.out, txt.file)
