# Script: 
# test gene level CN matrices for SCANB and BASIS for singificant differences in CN gain/loss alterations

#TODO: 

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

#-------------------------------------------------------------------------------
# packages
source("scripts/4_CN/src/cn_functions.R")
library(ggplot2)
library(tidyverse)
library(reshape2)
#library(data.table)
#library(purrr)
#library(readxl)

#-------------------------------------------------------------------------------
# set/create output directories
# for plots
output.path <- "output/plots/4_CN/"
dir.create(output.path)
# for data
data.path <- "data/SCANB/4_CN/processed/"
dir.create(data.path)
#-------------------------------------------------------------------------------
# input & output file paths
# input
scanb.gene.cna <- "data/SCANB/4_CN/processed/CNA_genelevel.RData"
basis.gene.cna <- "data/BASIS/4_CN/processed/CNA_genelevel_all.RData"
basis.anno <- "data/BASIS/1_clinical/raw/Summarized_Annotations_BASIS.RData"
# output
plot.file <- "output/plots/4_CN/COMBINED_HER2n_signifgenes.pdf"
txt.file <- "output/plots/4_CN/COMBINED_HER2n_signifgenes.txt"
signif.genes <- "data/COMBINED/4_CN/processed/CNA_genelevel.RData"
#-------------------------------------------------------------------------------
# storing objects 
plot.list <- list() # object to store plots; note: if the output is not in string format use capture.output()
txt.out <- c() # object to store text output
# snippet to save plots (paste to end of script)
# pdf(file = plot.file, onefile = TRUE, height = 5, width = 5)
# for(i in 1:length(plot.list)) { 
#   print(i)
#   print(plot.list[[i]])
# }
# dev.off()
# snippet to save text output (paste to end of script)
#writeLines(txt.out, txt.file)
#-------------------------------------------------------------------------------
# snipped to take exection time
#start.time <- Sys.time()
# CODE
#end.time <- Sys.time()
#time.taken <- end.time - start.time
#time.taken

################################################################################
# required data
################################################################################
scanb.cna <- loadRData(scanb.gene.cna)[["gainloss"]]
basis.cna <- loadRData(basis.gene.cna) %>% 
  dplyr::select(-c(chr,centerPos,Genome_pos))
#View(basis.cna)
basis.anno <- loadRData(basis.anno) %>% 
  filter(final.ER=="positive" & final.HER2=="negative") %>%
  dplyr::filter(PAM50_AIMS %in% c("LumA","LumB")) %>%
  dplyr::rename(PAM50=PAM50_AIMS) %>% 
  dplyr::select(sample_name,PAM50)

luma.ids <- basis.anno %>% filter(PAM50 == "LumA") %>% pull(sample_name)
lumb.ids <- basis.anno %>% filter(PAM50 == "LumB") %>% pull(sample_name)

#setdiff(luma.ids,names(basis.cna)) #1 luma sample doesnt have wgs data

# include only common genes
scanb.cna <- scanb.cna %>% drop_na() 
basis.cna <- basis.cna %>% drop_na()
common.genes <- intersect(scanb.cna$gene,basis.cna$gene)
scanb.cna <- scanb.cna %>% filter(gene %in% common.genes)
basis.cna <- basis.cna %>% filter(gene %in% common.genes)

# add gene annotation to basis from scanb
basis.cna$centerPos <- scanb.cna$centerPos[match(basis.cna$gene, scanb.cna$gene)]
basis.cna$Genome_pos <- scanb.cna$Genome_pos[match(basis.cna$gene, scanb.cna$gene)]
basis.cna$chr <- scanb.cna$chr[match(basis.cna$gene, scanb.cna$gene)]
basis.cna <- basis.cna %>% relocate(c(chr,centerPos,Genome_pos), .after=gene)

################################################################################

# save obj
# gene LumA.Gain.pval LumB.Gain.pval LumA.Loss.pval LumB.Loss.pval
LumA.Gain.pval <- c()
LumB.Gain.pval <- c()
LumA.Loss.pval <- c()
LumB.Loss.pval <- c()

pb = txtProgressBar(min = 0, max = length(common.genes), initial = 0, style = 3)
for(i in 1:length(common.genes)) {
  setTxtProgressBar(pb,i)
  # get data
  gene <- common.genes[i]
  #print(i)
  #print(gene)
  basis.gene.dat <- basis.cna[basis.cna$gene==gene,][,!names(basis.cna) %in% c(
    "gene","chr","centerPos","Genome_pos")]
  # 
  her2e.dat <- as.numeric(scanb.cna[scanb.cna$gene==gene,][,!names(scanb.cna) %in% c("gene","chr","centerPos","Genome_pos")])
  luma.dat <- as.numeric(basis.gene.dat[,names(basis.gene.dat) %in% luma.ids]) 
  lumb.dat <- as.numeric(basis.gene.dat[,names(basis.gene.dat) %in% lumb.ids]) 
  
  comp.list <- list("LumA"=luma.dat,"LumB"=lumb.dat)
  # test for all pairs
  for(j in 1:length(comp.list)) {
    comp.dat <- comp.list[[j]] 
    comp.group <- names(comp.list)[j]
    # get cross tables
    gain.tbl <- data.frame(her2e=c(sum(her2e.dat>0),sum(her2e.dat<1)), 
                           comp=c(sum(comp.dat>0),sum(comp.dat<1)), 
                           row.names = c("gain","no_gain"))
    loss.tbl <- data.frame(her2e=c(sum(her2e.dat<0),sum(her2e.dat>-1)), 
                           comp=c(sum(comp.dat<0),sum(comp.dat>-1)), 
                           row.names = c("loss","no_loss"))
    # test gain
    if(min(chisq.test(gain.tbl)$expected)<5) {
      gain.pval <- fisher.test(gain.tbl)$p.value
    } else {
      gain.pval <- chisq.test(gain.tbl)$p.value 
    }
    # test loss
    if(min(chisq.test(loss.tbl)$expected)<5) {
      loss.pval <- fisher.test(loss.tbl)$p.value
    } else {
      loss.pval <- chisq.test(loss.tbl)$p.value 
    }
    # save in vectors
    if (comp.group=="LumA") {
      LumA.Gain.pval[i] <- gain.pval
      LumA.Loss.pval[i] <- loss.pval
    } else if (comp.group=="LumB") {
      LumB.Gain.pval[i] <- gain.pval
      LumB.Loss.pval[i] <- loss.pval
    }
    
  }
  close(pb)
}


gene.test.df <- data.frame("gene"=common.genes, 
           "LumA.Gain.pval"=LumA.Gain.pval,
           "LumB.Gain.pval"=LumB.Gain.pval,
           "LumA.Loss.pval"=LumA.Loss.pval,
           "LumB.Loss.pval"=LumB.Loss.pval)

# adjust pval
gene.test.df$LumA.Gain.padj <- p.adjust(gene.test.df$LumA.Gain.pval, method = "fdr") # "bonferroni", "fdr"
gene.test.df$LumB.Gain.padj <- p.adjust(gene.test.df$LumB.Gain.pval, method = "fdr")
gene.test.df$LumA.Loss.padj <- p.adjust(gene.test.df$LumA.Loss.pval, method = "fdr")
gene.test.df$LumB.Loss.padj <- p.adjust(gene.test.df$LumB.Loss.pval, method = "fdr")
View(gene.test.df)
length(gene.test.df %>% filter(LumA.Gain.padj<=0.05) %>% pull(gene))
length(gene.test.df %>% filter(LumB.Gain.padj<=0.05) %>% pull(gene))
length(gene.test.df %>% filter(LumA.Loss.padj<=0.05) %>% pull(gene))
length(gene.test.df %>% filter(LumB.Loss.padj<=0.05) %>% pull(gene))
# how cn the same gene be signif lost and gained in luma lumb

save(gene.test.df, file= signif.genes)
