# Script: Test significnace of CNA freqs for each gene, create file
# Author: Lennart Hohmann
# Date: 20.02.2024
#-------------------
# empty environment 
rm(list=ls())
# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")
# cohort
cohort <- "SCANB"
#-------------------
# packages
source("./scripts/src/general_functions.R")
#source("./scripts/3_WGS/src/wgs_functions.R")
if (!require("pacman")) install.packages("pacman")
#pacman::p_load()
#-------------------
# set/create output directories
# for plots
output.path <- "./output/plots/4_CN/"
dir.create(output.path)
# for data
data.path <- "./data/SCANB/4_CN/processed/"
dir.create(data.path)
#-------------------
# input paths
infile.1 <- "./data/SCANB/4_CN/processed/CNA_genelevel_all.RData"
infile.2 <- "./data/SCANB/4_CN/processed/CNA_GLFreqs_all.RData"
# output paths
outfile.1 <- paste0(data.path,"CNA_genetest.RData")
#plot.file <- paste0(output.path,cohort,"_i.pdf")
#txt.file <- paste0(output.path,cohort,"_i.txt")
#-------------------
# storing objects 
#plot.list <- list() # object to store plots
#plot.parameters <- list() # object to store parameters to plot base R plots again later
#txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

ascat.df.all <- loadRData(infile.1)[["ascat.df.all"]]
ascat.df.all[,2:ncol(ascat.df.all)] <- lapply(ascat.df.all[,2:ncol(ascat.df.all)], as.numeric)
her2e.ids <- loadRData(infile.1)[["subtype.samples"]][["her2e"]]
luma.ids <- loadRData(infile.1)[["subtype.samples"]][["LumA"]]
lumb.ids <- loadRData(infile.1)[["subtype.samples"]][["LumB"]]

#######################################################################
# test genes Loss and Gain freqs
#######################################################################

LumA.Loss.pval <- c()
LumA.Gain.pval <- c()
LumB.Loss.pval <- c()
LumB.Gain.pval <- c()

pb = txtProgressBar(min = 0, max = length(ascat.df.all$gene), initial = 0, style = 3)
for(i in 1:length(ascat.df.all$gene)) {
  setTxtProgressBar(pb,i)
  # gene metadata
  gene <- ascat.df.all$gene[i]
  gene.anno <- ascat.df.all[ascat.df.all$gene==gene,
                            c("gene","chr","start","end")]
  
  # # gene CNA data with PAM50 sample annotation
  # gene.dat <- as.data.frame(
  #   t(ascat.df.all[ascat.df.all$gene==gene,
  #                  !names(ascat.df.all) %in% c("gene","chr","start","end")]))
  # gene.dat$sampleID <- rownames(gene.dat)
  # gene.dat$PAM50 <- ifelse(gene.dat$sampleID %in% her2e.ids, "her2e",
  #                         ifelse(gene.dat$sampleID %in% luma.ids, "LumA",
  #                                ifelse(gene.dat$sampleID %in% lumb.ids, "LumB", NA)))
  # gene.dat$sampleID <- NULL
  # names(gene.dat)[1] <- "gene.CNA"
  
  #View(gene.dat)
  her2e.dat <- unlist(ascat.df.all[ascat.df.all$gene==gene,her2e.ids])
  luma.dat <- unlist(ascat.df.all[ascat.df.all$gene==gene,luma.ids])
  lumb.dat <- unlist(ascat.df.all[ascat.df.all$gene==gene,lumb.ids])

  comp.list <- list("LumA"=luma.dat,"LumB"=lumb.dat)
  
  for(j in 1:length(comp.list)) {
    comp.dat <- comp.list[[j]] 
    comp.group <- names(comp.list)[j]
    
    gain.tbl <- data.frame(her2e=c(sum(her2e.dat==1),sum(her2e.dat!=1)), 
                           comp=c(sum(comp.dat==1),sum(comp.dat!=1)), 
                           row.names = c("gain","no_gain"))
    
    loss.tbl <- data.frame(her2e=c(sum(her2e.dat==-1),sum(her2e.dat!=-1)), 
                           comp=c(sum(comp.dat==-1),sum(comp.dat!=-1)), 
                           row.names = c("loss","no_loss"))
    
    # test gain
    gain.pval <- fisher.test(gain.tbl)$p.value
    # test loss
    loss.pval <- fisher.test(loss.tbl)$p.value
    
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



gene.test.df <- data.frame("gene"=ascat.df.all$gene, 
                           "LumA.Gain.pval"=LumA.Gain.pval,
                           "LumB.Gain.pval"=LumB.Gain.pval,
                           "LumA.Loss.pval"=LumA.Loss.pval,
                           "LumB.Loss.pval"=LumB.Loss.pval)


#save(gene.test.df, file= signif.genes)

# adjust pval
gene.test.df$LumA.Gain.padj <- p.adjust(gene.test.df$LumA.Gain.pval, method = "fdr") #fdr
gene.test.df$LumB.Gain.padj <- p.adjust(gene.test.df$LumB.Gain.pval, method = "fdr")
gene.test.df$LumA.Loss.padj <- p.adjust(gene.test.df$LumA.Loss.pval, method = "fdr")
gene.test.df$LumB.Loss.padj <- p.adjust(gene.test.df$LumB.Loss.pval, method = "fdr")
#View(gene.test.df)
#length(gene.test.df[gene.test.df$LumA.Gain.padj<=0.05,]$gene)
#length(gene.test.df[gene.test.df$LumB.Gain.padj<=0.05,]$gene)
#length(gene.test.df[gene.test.df$LumA.Loss.padj<=0.05,]$gene)
#length(gene.test.df[gene.test.df$LumB.Loss.padj<=0.05,]$gene)

save(gene.test.df, file= outfile.1)
