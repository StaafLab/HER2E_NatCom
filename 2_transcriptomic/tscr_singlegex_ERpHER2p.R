# Script: In HER2-positive: Single gene expression assessment (incl. ERBB2, ESR1, FGFR4)

#TODO

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

# output filenames
plot.list <- list() # object to store plots; note: if the output is not in string format use capture.output()
plot.file <- paste(output.path,cohort,"_HER2p_singlegex_plots.pdf", sep="")

txt.out <- c() # object to store text output
txt.file <- paste(output.path,cohort,"_HER2p_singlegex_text.txt", sep="")

#packages
source("./scripts/2_transcriptomic/src/tscr_functions.R")
library(ggplot2)
library(tidyverse)
library(readxl)
library(biomaRt)
library(gridExtra)
library(ggsignif)
library(janitor)
library(scales)

#######################################################################
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
             c("HER2n_HER2E","HER2p_HER2E","HER2p_nonHER2E"))
  
  # load gex data
  gex.data <- scanb_gex_load(gex.path = "data/SCANB/2_transcriptomic/raw/genematrix_noNeg.Rdata", geneanno.path = "data/SCANB/1_clinical/raw/Gene.ID.ann.Rdata", ID.type = "Gene.Name") %>%
    dplyr::select(any_of(anno$sampleID)) %>% # select subgroup gex
    select_if(~ !any(is.na(.))) # otherwise error when scaling 
  
  # log transform FPKM data
  gex.data <- as.data.frame(log2(gex.data + 1))
  # z-transform
  gex.data <- as.data.frame(t(apply(gex.data, 1, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y))))) # for some rows there may be 0 variance so i have to handle these cases
  
  # source file export
  source.dat <- as.data.frame(t(gex.data[c("ESR1","ERBB2","FGFR4"),]))
  source.dat$sampleID <- rownames(source.dat)
  rownames(source.dat) <- NULL
  source.dat <- source.dat[c("sampleID","ESR1","ERBB2","FGFR4")]
  source.dat <- merge(anno[c("sampleID","Group")],source.dat,by="sampleID")
  save(source.dat,file = "./output/source_data/R_objects/Figure_5_singlegexHER2p.RData")
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
# look into the expression of selected genes # CONTINUE HERE
#######################################################################
# *<0.05, **<0.01 ***<0.001 ****< 0.0001 ns not significant

gene.vec <- c("ERBB2","ESR1","FGFR4")

for (gene in gene.vec) {

  gene.gex <- get_gex_hp(gene,gex.data,anno)
  
  # base statistics
  stats <- capture.output(get_stats(gene.gex,"Group",gene))
  
  # comparison data
  HER2n_HER2E.dat <- gene.gex[gene.gex$Group=="HER2n_HER2E",][[gene]] # unlist()
  HER2p_HER2E.dat <- gene.gex[gene.gex$Group=="HER2p_HER2E",][[gene]]
  HER2p_nonHER2E.dat <- gene.gex[gene.gex$Group=="HER2p_nonHER2E",][[gene]]
  
  # pair comp 1
  res <- mwu_test(HER2n_HER2E.dat,HER2p_HER2E.dat)
  txt.out <- append(txt.out,
                    c(gene," statistics: HER2n_HER2E.dat vs. HER2p_HER2E.dat",
                      capture.output(res)))
  
  # pair comp 2
  res <- mwu_test(HER2n_HER2E.dat,HER2p_nonHER2E.dat)
  txt.out <- append(txt.out,
                    c(gene," statistics: HER2n_HER2E.dat vs. HER2p_nonHER2E.dat",
                      capture.output(res)))
  
  # plot
  plot <- three_boxplot(gene.gex,
                        group.var = "Group",
                        test.var = gene,
                        g1="HER2n_HER2E",g2="HER2p_nonHER2E",g3="HER2p_HER2E",
                        colors=setNames(c("#d334eb","#d8b365","#5ab4ac"),
                                        c("HER2n_HER2E","HER2p_nonHER2E", 
                                          "HER2p_HER2E")),
                        ylab = "Expression (log2)", 
                        #ylim = if (cohort=="SCANB") {c(-6.5,9)} else {c(-3,3)},
                        title = paste(gene," expression in HER2/PAM50 subtypes", 
                        sep=""))
  
  plot.list <- append(plot.list, list(plot))
}
#######################################################################

# save plots
pdf(file = plot.file, 
    onefile = TRUE, width = 15, height = 15) 
for (i in 1:length(plot.list)) {
  print(plot.list[[i]])
}
dev.off()

# save text output
writeLines(txt.out, txt.file)
