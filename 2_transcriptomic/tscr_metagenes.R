# Script: Metagene analysis in SCAN-B and METABRIC

# TODO:
# try once with scaling MB data and without and compare output

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
plot.file <- paste(output.path,cohort,"_HER2n_metagenes_plots.pdf", sep="")
txt.out <- c() # object to store text output
txt.file <- paste(output.path,cohort,"_HER2n_metagenes_text.txt", sep="")

#packages
source("scripts/2_transcriptomic/src/tscr_functions.R")
library(ggplot2)
library(tidyverse)
#library(matrixStats)
library(readxl)
library(biomaRt)
library(gridExtra)
library(ggsignif)
library(janitor)
library(data.table)

#######################################################################
# load metagene definitions
#######################################################################

# metagene definitions
metagene.def <- as.data.frame(read_excel(
    "data/SCANB/2_transcriptomic/raw/metagene_definitions.XLSX")) %>% 
    dplyr::rename(entrezgene_id = `Entrez Gene ID`,
                  module = `Module Name`,
                  gene_symbol = `Gene symbol`)

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
    gex.data <- scanb_gex_load(gex.path = "data/SCANB/2_transcriptomic/raw/genematrix_noNeg.Rdata", geneanno.path = "data/SCANB/1_clinical/raw/Gene.ID.ann.Rdata", ID.type = "EntrezGene") %>%
      dplyr::select(any_of(anno$sampleID)) %>% # select subgroup gex
      filter(row.names(.) %in% metagene.def$entrezgene_id) %>% 
      select_if(~ !any(is.na(.))) # otherwise error when scaling 

    # log transform FPKM data
    gex.data <- as.data.frame(log2(gex.data + 1))
    # z-transform
    gex.data <- as.data.frame(t(apply(gex.data, 1, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y))))) # for some rows there may be 0 variance so i have to handle these cases
    
#-----------------------------------------------------------------------#
    
} else if (cohort=="METABRIC") {
    
  #load annotation data
  load("data/METABRIC/1_clinical/raw/Merged_annotations.RData")
  
  # extract relevant variables
  anno <- loadRData("data/METABRIC/1_clinical/processed/Merged_annotations_ERpHER2n.RData") %>% 
    dplyr::rename(sampleID=METABRIC_ID) # rename to match SCANB variables
  
  # load and select subgroup data
  gex.data <- metabric_gex_load("./data/METABRIC/2_transcriptomic/raw/data_mRNA_median_all_sample_Zscores.txt",ID.type = "Entrez_Gene_Id") %>% 
    dplyr::select(any_of(anno$sampleID)) %>% 
    mutate_all(function(x) as.numeric(x)) %>% 
    filter(row.names(.) %in% metagene.def$entrezgene_id) %>% 
    select_if(~ !any(is.na(.))) 
    
    # exclude samples from anno without associated gex data
    anno <- anno %>% 
        filter(sampleID %in% colnames(gex.data))
}


# check if genes were not in the scanb data
setdiff(metagene.def$entrezgene_id,rownames(gex.data)) 

#METABRIC: 125 -> ; 9370 -> ; 1520 -> ; 3507 -> "IGHM" ; 28755 <- "TRAC"
#SCANB: 3507 -> "IGHM" ; 28755 <- "TRAC"

#######################################################################
# 4. calc. score for each metagene in each sample
#######################################################################

#basal
basal.scores <- mgscore(metagene = "Basal",
                metagene.def = metagene.def,
                gex.data = gex.data) 

#earlyresponse
earlyresponse.scores <- mgscore(metagene = "Early_response",
                        metagene.def = metagene.def,
                        gex.data = gex.data)

metagene.scores <- merge(basal.scores,earlyresponse.scores,by=0) %>% 
    column_to_rownames(var = "Row.names")

#immuneresponse
ir.scores <- mgscore(metagene = "IR",
                                metagene.def = metagene.def,
                                gex.data = gex.data)

metagene.scores <- merge(metagene.scores,ir.scores,by=0) %>% column_to_rownames(var = "Row.names")

#lipid
lipid.scores <- mgscore(metagene = "Lipid",
                        metagene.def = metagene.def,
                        gex.data = gex.data)

metagene.scores <- merge(metagene.scores,lipid.scores,by=0) %>% column_to_rownames(var = "Row.names")

#mitoticcheckpoint
mitoticcheckpoint.scores <- mgscore(metagene = "Mitotic_checkpoint",
                                    metagene.def = metagene.def,
                                    gex.data = gex.data)

metagene.scores <- merge(metagene.scores,mitoticcheckpoint.scores,by=0) %>% column_to_rownames(var = "Row.names")


#mitoticprogression
mitoticprogression.scores <- mgscore(metagene = "Mitotic_progression",
                                    metagene.def = metagene.def,
                                    gex.data = gex.data)

metagene.scores <- merge(metagene.scores,mitoticprogression.scores,by=0) %>% column_to_rownames(var = "Row.names")


#SR
sr.scores <- mgscore(metagene = "SR",
                    metagene.def = metagene.def,
                    gex.data = gex.data)

metagene.scores <- merge(metagene.scores,sr.scores,by=0) %>% column_to_rownames(var = "Row.names")

#stroma
stroma.scores <- mgscore(metagene = "Stroma",
                     metagene.def = metagene.def,
                     gex.data = gex.data)

metagene.scores <- merge(metagene.scores,stroma.scores,by=0) %>% column_to_rownames(var = "Row.names")

#######################################################################
# 5. statistics for metagene scores between groups
#######################################################################

mg.pvals <- data.frame()
metagene.vec <- colnames(metagene.scores)
for (metagene in metagene.vec) {
  
  metagene.score <- metagene.scores[metagene]
  
  # comparison data
  Her2.dat <- metagene.score[anno[anno$PAM50=="Her2",]$sampleID,]
  Lumb.dat <- metagene.score[anno[anno$PAM50=="LumB",]$sampleID,]
  LumA.dat <- metagene.score[anno[anno$PAM50=="LumA",]$sampleID,]
  
  # pair comp 1
  res <- mwu_test(Her2.dat,Lumb.dat)
  Lumb.pval <- res$p.value
  txt.out <- append(txt.out,
                    c(metagene," statistics: Her2.dat vs. Lumb.dat",
                      capture.output(res)))
  
  # pair comp 2
  res <- mwu_test(Her2.dat,LumA.dat)
  LumA.pval <- res$p.value
  txt.out <- append(txt.out,
                    c(metagene," statistics: Her2.dat vs. LumA.dat",
                      capture.output(res)))
  
  # append to final res df
  mg.pvals[
    metagene,
    c("Her2.Lumb.pval","Her2.LumA.pval")] <- c(Lumb.pval,LumA.pval)
}

# make combined data and anno object for plotting
mg.anno <- merge(metagene.scores %>% rownames_to_column(var="sampleID"),anno[,c("sampleID","PAM50")],by="sampleID")

mg.anno.list <- list(mg.anno, mg.pvals)
save(mg.anno.list,file = paste(data.path,"mg_anno.RData",sep=""))

mg.pvals$Her2.Lumb.padj <- p.adjust(mg.pvals$Her2.Lumb.pval, 
                                    method = p.adjust.methods, 
                                    n = 16)
mg.pvals$Her2.Luma.padj <- p.adjust(mg.pvals$Her2.LumA.pval, 
                                    method = p.adjust.methods, 
                                    n = 16)


mg.pvals <- mg.pvals %>% 
  mutate_if(is.numeric, round, digits=4) %>% 
  rownames_to_column(var="Metagene")
txt.out <- append(txt.out, c("FDR ADJUSTED P-VALUES"))
txt.out <- append(txt.out, c(capture.output(mg.pvals)))

#######################################################################
# 5. Boxplots
#######################################################################
# *<0.05, **<0.01 ***<0.001 ****< 0.0001   ns not significant

# SCANB
plot.list <- append(plot.list, list(
    three_boxplot(mg.anno,
              group.var = "PAM50",
              test.var = "Basal",
              g1="Her2",g2="LumA",g3="LumB",
              colors=setNames(c("#d334eb","#2176d5","#34c6eb"),
                              c("Her2","LumA","LumB")),
              ylim = if (cohort=="SCANB") {c(-1.5,3.5)} else {c(-1.5,2.5)},
              ylab = "Metagene score", 
              title = "Basal metagene scores in PAM50 subtypes (ERpHER2n)")))

plot.list <- append(plot.list, list(
    three_boxplot(mg.anno,
                  group.var = "PAM50",
                  test.var = "Early_response",
                  g1="Her2",g2="LumA",g3="LumB",
                  colors=setNames(c("#d334eb","#2176d5","#34c6eb"),
                                  c("Her2","LumA","LumB")),
                  ylim = if (cohort=="SCANB") {c(-2.5,4.2)} else {c(-2.5,3)},
                  ylab = "Metagene score", 
                  title = "Early_response metagene scores in PAM50 subtypes (ERpHER2n)")))

plot.list <- append(plot.list, list(
    three_boxplot(mg.anno,
                  group.var = "PAM50",
                  test.var = "IR",
                  g1="Her2",g2="LumA",g3="LumB",
                  colors=setNames(c("#d334eb","#2176d5","#34c6eb"),
                                  c("Her2","LumA","LumB")),
                  ylim = if (cohort=="SCANB") {c(-2,5)} else {c(-2,3.5)},
                  ylab = "Metagene score", 
                  title = "IR metagene scores in PAM50 subtypes (ERpHER2n)")))
  
plot.list <- append(plot.list, list(
    three_boxplot(mg.anno,
                  group.var = "PAM50",
                  test.var = "Lipid",
                  g1="Her2",g2="LumA",g3="LumB",
                  colors=setNames(c("#d334eb","#2176d5","#34c6eb"),
                                  c("Her2","LumA","LumB")),
                  ylim = if (cohort=="SCANB") {c(-2,3.5)} else {c(-1.5,3)},
                  ylab = "Metagene score", 
                  title = "Lipid metagene scores in PAM50 subtypes (ERpHER2n)")))

plot.list <- append(plot.list, list(
    three_boxplot(mg.anno,
                  group.var = "PAM50",
                  test.var = "Mitotic_checkpoint",
                  g1="Her2",g2="LumA",g3="LumB",
                  colors=setNames(c("#d334eb","#2176d5","#34c6eb"),
                                  c("Her2","LumA","LumB")),
                  ylim = if (cohort=="SCANB") {c(-2,4.5)} else {c(-2,2.5)},
                  ylab = "Metagene score", 
                  title = "Mitotic_checkpoint metagene scores in PAM50 subtypes (ERpHER2n)")))
    
plot.list <- append(plot.list, list(
    three_boxplot(mg.anno,
                  group.var = "PAM50",
                  test.var = "Mitotic_progression",
                  g1="Her2",g2="LumA",g3="LumB",
                  colors=setNames(c("#d334eb","#2176d5","#34c6eb"),
                                  c("Her2","LumA","LumB")),
                  ylim = if (cohort=="SCANB") {c(-2,4.5)} else {c(-2,3)},
                  ylab = "Metagene score", 
                  title = "Mitotic_progression metagene scores in PAM50 subtypes (ERpHER2n)")))
    
plot.list <- append(plot.list, list(
    three_boxplot(mg.anno,
                  group.var = "PAM50",
                  test.var = "SR",
                  g1="Her2",g2="LumA",g3="LumB",
                  colors=setNames(c("#d334eb","#2176d5","#34c6eb"),
                                  c("Her2","LumA","LumB")),
                  ylim = if (cohort=="SCANB") {c(-2.5,2.5)} else {c(-1.5,2)},
                  ylab = "Metagene score", 
                  title = "SR metagene scores in PAM50 subtypes (ERpHER2n)")))
  
plot.list <- append(plot.list, list(
    three_boxplot(mg.anno,
                  group.var = "PAM50",
                  test.var = "Stroma",
                  g1="Her2",g2="LumA",g3="LumB",
                  colors=setNames(c("#d334eb","#2176d5","#34c6eb"),
                                  c("Her2","LumA","LumB")),
                  ylim = if (cohort=="SCANB") {c(-3.5,3)} else {c(-3,3)},
                  ylab = "Metagene score", 
                  title = "Stroma metagene scores in PAM50 subtypes (ERpHER2n)")))

############################################################################

# save plots
pdf(file = plot.file, 
    onefile = TRUE, width = 15, height = 15) 
for (i in 1:length(plot.list)) {
  print(plot.list[[i]])
}
dev.off()

# save text output
writeLines(txt.out, txt.file)
