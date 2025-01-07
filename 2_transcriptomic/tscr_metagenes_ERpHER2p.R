# Script: Metagene analysis for ERpHER2nHER2E and ERpHER2p contrast groups in SCAN-B and METABRIC

# TODO:

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

txt.out <- c() # object to store text output
txt.file <- paste(output.path,cohort,"_HER2p_metagenes_text.txt", sep="")


#packages
source("scripts/2_transcriptomic/src/tscr_functions.R")
library(ggplot2)
library(tidyverse)
#library(matrixStats)
library(readxl)
library(biomaRt)
library(gridExtra)
library(ggsignif)
library(forcats)
library(janitor)

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
    gex.data <- metabric_gex_load("./data/METABRIC/2_transcriptomic/raw/data_mRNA_median_all_sample_Zscores.txt",ID.type = "Entrez_Gene_Id") %>% 
      dplyr::select(any_of(anno$sampleID)) %>% 
      mutate_all(function(x) as.numeric(x)) %>% 
      filter(row.names(.) %in% metagene.def$entrezgene_id) %>% 
      select_if(~ !any(is.na(.))) 
    
    # exclude samples from anno without associated gex data
    anno <- anno %>% 
        filter(sampleID %in% colnames(gex.data))
    
}
#View(head(gex.data))
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
  HER2n_HER2E.dat <- metagene.score[anno[anno$Group=="HER2n_HER2E",]$sampleID,]
  HER2p_HER2E.dat <- metagene.score[anno[anno$Group=="HER2p_HER2E",]$sampleID,]
  HER2p_nonHER2E.dat <- metagene.score[anno[anno$Group=="HER2p_nonHER2E",]$sampleID,]
  
  # pair comp 1
  res <- mwu_test(HER2n_HER2E.dat,HER2p_HER2E.dat)
  HER2p_HER2E.pval <- res$p.value
  txt.out <- append(txt.out,
                    c(metagene," statistics: HER2n_HER2E.dat vs. HER2p_HER2E.dat",
                      capture.output(res)))
  
  # pair comp 2
  res <- mwu_test(HER2n_HER2E.dat,HER2p_nonHER2E.dat)
  HER2p_nonHER2E.pval <- res$p.value
  txt.out <- append(txt.out,
                    c(metagene," statistics: HER2n_HER2E.dat vs. HER2p_nonHER2E.dat",
                      capture.output(res)))
  
  # append to final res df
  mg.pvals[
    metagene,
    c("HER2n_HER2E.HER2p_HER2E.pval","HER2n_HER2E.HER2p_nonHER2E.pval")] <- c(HER2p_HER2E.pval,HER2p_nonHER2E.pval)
}

# make combined data and anno object for plotting
mg.anno <- merge(metagene.scores %>% 
                   rownames_to_column(var="sampleID"),anno[,c("sampleID","Group")],by="sampleID")

mg.anno.list <- list(mg.anno, mg.pvals)
save(mg.anno.list,file = paste(data.path,"mg_anno_HER2p.RData",sep=""))


mg.pvals$HER2n_HER2E.HER2p_HER2E.padj <- p.adjust(mg.pvals$HER2n_HER2E.HER2p_HER2E.pval, 
                                    method = p.adjust.methods, 
                                    n = 16)
mg.pvals$HER2n_HER2E.HER2p_nonHER2E.padj <- p.adjust(mg.pvals$HER2n_HER2E.HER2p_nonHER2E.pval, 
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
# plot
# *<0.05, **<0.01 ***<0.001 ****< 0.0001   ns not significant

# list to save plots
plot.list <- list()


plot.list <- append(plot.list, list(
  three_boxplot(mg.anno,
                group.var = "Group",
                test.var = "Basal",
                g1="HER2n_HER2E",g2="HER2p_nonHER2E",g3="HER2p_HER2E",
                colors=setNames(c("#d334eb","#d8b365","#5ab4ac"),
                                c("HER2n_HER2E","HER2p_nonHER2E","HER2p_HER2E")),
                ylim = if (cohort=="SCANB") {c(-1.5,4)} else {c(-1.5,4)},
                ylab = "Metagene score", 
                title = "Basal metagene scores in HER2/PAM50 subtypes")))

plot.list <- append(plot.list, list(
  three_boxplot(mg.anno,
                group.var = "Group",
                test.var = "IR",
                g1="HER2n_HER2E",g2="HER2p_nonHER2E",g3="HER2p_HER2E",
                colors=setNames(c("#d334eb","#d8b365","#5ab4ac"),
                                c("HER2n_HER2E","HER2p_nonHER2E","HER2p_HER2E")),
                ylim = if (cohort=="SCANB") {c(-2,4.1)} else {c(-2,4.1)},
                ylab = "Metagene score", 
                title = "IR metagene scores in HER2/PAM50 subtypes")))
    
plot.list <- append(plot.list, list(
  three_boxplot(mg.anno,
                group.var = "Group",
                test.var = "Lipid",
                g1="HER2n_HER2E",g2="HER2p_nonHER2E",g3="HER2p_HER2E",
                colors=setNames(c("#d334eb","#d8b365","#5ab4ac"),
                                c("HER2n_HER2E","HER2p_nonHER2E","HER2p_HER2E")),
                ylim = if (cohort=="SCANB") {c(-1.8,3.5)} else {c(-1.8,3.5)},
                ylab = "Metagene score", 
                title = "Lipid metagene scores in HER2/PAM50 subtypes")))
                                        
plot.list <- append(plot.list, list(
  three_boxplot(mg.anno,
                group.var = "Group",
                test.var = "Mitotic_checkpoint",
                g1="HER2n_HER2E",g2="HER2p_nonHER2E",g3="HER2p_HER2E",
                colors=setNames(c("#d334eb","#d8b365","#5ab4ac"),
                                c("HER2n_HER2E","HER2p_nonHER2E","HER2p_HER2E")),
                ylim = if (cohort=="SCANB") {c(-3,3.8)} else {c(-3,3.8)},
                ylab = "Metagene score", 
                title = "Mitotic_checkpoint metagene scores in HER2/PAM50 subtypes")))

plot.list <- append(plot.list, list(
  three_boxplot(mg.anno,
                group.var = "Group",
                test.var = "Mitotic_progression",
                g1="HER2n_HER2E",g2="HER2p_nonHER2E",g3="HER2p_HER2E",
                colors=setNames(c("#d334eb","#d8b365","#5ab4ac"),
                                c("HER2n_HER2E","HER2p_nonHER2E","HER2p_HER2E")),
                ylim = if (cohort=="SCANB") {c(-3,4)} else {c(-3,4)},
                ylab = "Metagene score", 
                title = "Mitotic_progression metagene scores in HER2/PAM50 subtypes")))


plot.list <- append(plot.list, list(
  three_boxplot(mg.anno,
                group.var = "Group",
                test.var = "SR",
                g1="HER2n_HER2E",g2="HER2p_nonHER2E",g3="HER2p_HER2E",
                colors=setNames(c("#d334eb","#d8b365","#5ab4ac"),
                                c("HER2n_HER2E","HER2p_nonHER2E","HER2p_HER2E")),
                ylim = if (cohort=="SCANB") {c(-4,3)} else {c(-4,3)},
                ylab = "Metagene score", 
                title = "SR metagene scores in HER2/PAM50 subtypes")))

plot.list <- append(plot.list, list(
  three_boxplot(mg.anno,
                group.var = "Group",
                test.var = "Stroma",
                g1="HER2n_HER2E",g2="HER2p_nonHER2E",g3="HER2p_HER2E",
                colors=setNames(c("#d334eb","#d8b365","#5ab4ac"),
                                c("HER2n_HER2E","HER2p_nonHER2E","HER2p_HER2E")),
                ylim = if (cohort=="SCANB") {c(-2.7,2.7)} else {c(-2.7,2.7)},
                ylab = "Metagene score", 
                title = "Stroma metagene scores in HER2/PAM50 subtypes")))

###########################################################################

#plot
pdf(file = paste(output.path,cohort,"_HER2p_metagenes.pdf", sep=""), 
    onefile = TRUE, width = 15, height = 15) 

for (i in 1:length(plot.list)) {
    print(plot.list[[i]])
}

dev.off()

# save text output
writeLines(txt.out, txt.file)

###########################################################################
