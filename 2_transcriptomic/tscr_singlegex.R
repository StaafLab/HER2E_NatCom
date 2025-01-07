# Script: Single gene expression assessment (incl. ERBB2, ESR1, FGFR4)

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
plot.file <- paste(output.path,cohort,"_HER2n_singlegex_plots.pdf", sep="")
txt.out <- c() # object to store text output
txt.file <- paste(output.path,cohort,"_HER2n_singlegex_text.txt", sep="")

# packages
source("./scripts/2_transcriptomic/src/tscr_functions.R")
library(ggplot2)
library(tidyverse)
#library(matrixStats)
library(readxl)
library(biomaRt)
library(gridExtra)
library(ggsignif)
library(janitor)
library(rstatix)
library(scales)

#######################################################################
# Cohort-specific data preprocessing including selection of  
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
  
  # source file export
  source.dat <- as.data.frame(t(gex.data[c("CD274","ESR1","ERBB2","FGFR4"),]))
  source.dat$sampleID <- rownames(source.dat)
  rownames(source.dat) <- NULL
  source.dat <- source.dat[c("sampleID","CD274","ESR1","ERBB2","FGFR4")]
  source.dat <- merge(anno[c("sampleID","PAM50")],source.dat,by="sampleID")
  save(source.dat,file = "./output/source_data/R_objects/Figure_3_singlegex.RData")
  #---------------------------------------------------------------------#
  
} else if (cohort=="METABRIC") {
  
  # load anno data
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
# look into the expression of selected genes
#######################################################################
# extreme <- gene.gex[is_outlier(gene.gex$ERBB2, coef = 3),]
#extreme %>% filter(ERBB2 > 0)
#c("S004036","S004794","S005551","S005605","S006395","S008332")

gene.vec <- c("ERBB2","ESR1","FGFR4","CD274") #TACSTD2

for (gene in gene.vec) {
  
  gene.gex <- get_gex(gene,gex.data,anno)
  
  # base statistics
  stats <- capture.output(get_stats(gene.gex,"PAM50",gene))
  txt.out <- append(txt.out,
                    c(gene,stats))
  # comparison data
  Her2.dat <- gene.gex[gene.gex$PAM50=="Her2",][[gene]] # unlist()
  LumB.dat <- gene.gex[gene.gex$PAM50=="LumB",][[gene]]
  LumA.dat <- gene.gex[gene.gex$PAM50=="LumA",][[gene]]
  
  # pair comp 1
  res <- mwu_test(Her2.dat,LumB.dat)
  txt.out <- append(txt.out,
                    c(gene," statistics: Her2.dat vs. LumB.dat",
                      capture.output(res)))
  
  # pair comp 2
  res <- mwu_test(Her2.dat,LumA.dat)
  txt.out <- append(txt.out,
                    c(gene," statistics: Her2.dat vs. LumA.dat",
                      capture.output(res)))
  
  # plot
  plot <- three_boxplot(gene.gex,
                        group.var = "PAM50",
                        test.var = gene,
                        g1="Her2",g2="LumA",g3="LumB",
                        colors=setNames(c("#d334eb","#2176d5","#34c6eb"),
                                        c("Her2","LumA","LumB")),
                        ylab = "Expression (log2)", 
                        title = paste(
                          gene," expression in PAM50 subtypes (ERpHER2n)", 
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
