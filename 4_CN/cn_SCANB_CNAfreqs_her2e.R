# Script: Create summarized CNA files for SCANB/BASIS and Calc. CNA Frequencies in SCAN-B
# Author: Lennart Hohmann
# Date: 20.02.2024
#-------------------
# empty environment
rm(list=ls())
# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")
# cohort
#cohort <- "SCANB"
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
#infile.1 <- "./data/SCANB/0_GroupSamples/ERpHER2n_sampleIDs.RData"
infile.2 <- "./data/SCANB/3_genomic/processed/ASCAT_genelevel.RData"
infile.7 <- "./data/BASIS/3_genomic/processed/ASCAT_genelevel.RData"
infile.8 <- "./data/BASIS/1_clinical/raw/Summarized_Annotations_BASIS.RData"
infile.9 <- "./data/SCANB/3_genomic/raw/HER2_enriched_June23_ForJohan.xlsx"
#infile.10 <- "./data/SCANB/3_WGS/raw/MergedAnnotations_ERp_Cohort_FailFiltered.RData"
# output paths
outfile.1 <- paste0(data.path,"CNA_genelevel_all.RData")
outfile.2 <- paste0(data.path,"CNA_GLFreqs_all.RData")
#plot.file <- paste0(output.path,cohort,"_i.pdf")
#txt.file <- paste0(output.path,cohort,"_i.txt")
#-------------------
# storing objects 
#plot.list <- list() # object to store plots
#plot.parameters <- list() # object to store parameters to plot base R plots again later
#txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load SCANB data - her2e subtype
#######################################################################

# wgs QC
qc.scanb <- as.data.frame(read_excel(infile.9, sheet = "Samples"))
pass.scanb.samples <- qc.scanb$Sample # all pass

# load her2e ids
her2e.ids <- pass.scanb.samples

# load ASCAT gene data
ascat.list.scanb <- loadRData(infile.2)
names(ascat.list.scanb) <- gsub("\\..*", "", names(ascat.list.scanb))
ascat.list.her2e <- ascat.list.scanb[names(ascat.list.scanb) %in% her2e.ids]
ascat.list.her2e <- ascat.list.scanb[pass.scanb.samples]

# get data format: gene sample1 sample2 ... 
ascat.df.her2e <- do.call(rbind, lapply(ascat.list.her2e, function(x) x$CNA))
ascat.df.her2e <- t(ascat.df.her2e)
ascat.df.her2e <- cbind(ascat.list.her2e[[1]][c("gene","chr","start","end")],ascat.df.her2e)
#View(ascat.df.her2e)

# get freqs
CNA.freqs.her2e <- ascat.df.her2e
# calc loss/gain freqs per group
CNA.freqs.her2e$freqloss.her2e <- apply(
  CNA.freqs.her2e[,5:ncol(CNA.freqs.her2e)], 1, function(x) (
    length(which(x==-1))/ncol(CNA.freqs.her2e[,5:ncol(CNA.freqs.her2e)]))*-100) # i add a minus to make it easier for plotting

CNA.freqs.her2e$freqgain.her2e <- apply(
  CNA.freqs.her2e[,5:ncol(CNA.freqs.her2e)], 1, function(x) (
    length(which(x==1))/ncol(CNA.freqs.her2e[,5:ncol(CNA.freqs.her2e)]))*100)

CNA.freqs.her2e <- CNA.freqs.her2e[c("gene","chr","start","end","freqloss.her2e","freqgain.her2e")]

#######################################################################
# load BASIS data - LumA and LumB subtypes
#######################################################################

# get relevant sample IDs
basis.anno <- loadRData(infile.8)
basis.anno <- basis.anno[basis.anno$ClinicalGroup == "ERposHER2neg" & 
                           basis.anno$PAM50_AIMS %in% c("LumA","LumB"),
                         c("sample_name","PAM50_AIMS")]

# LumA and B data
luma.ids <- basis.anno$sample_name[basis.anno$PAM50_AIMS=="LumA"]
lumb.ids <- basis.anno$sample_name[basis.anno$PAM50_AIMS=="LumB"]
ascat.list.basis <- loadRData(infile.7)
ascat.list.luma <- ascat.list.basis[names(ascat.list.basis) %in% luma.ids]
ascat.list.lumb <- ascat.list.basis[names(ascat.list.basis) %in% lumb.ids]

# prep luma files
ascat.df.luma <- do.call(rbind, lapply(ascat.list.luma, function(x) x$CNA))
ascat.df.luma <- t(ascat.df.luma)
ascat.df.luma <- cbind(ascat.list.luma[[1]][c("gene","chr","start","end")],ascat.df.luma)
# get freqs
CNA.freqs.luma <- ascat.df.luma
CNA.freqs.luma$freqloss.LumA <- apply(
  CNA.freqs.luma[,5:ncol(CNA.freqs.luma)], 1, function(x) (
    length(which(x==-1))/ncol(CNA.freqs.luma[,5:ncol(CNA.freqs.luma)]))*-100)
CNA.freqs.luma$freqgain.LumA <- apply(
  CNA.freqs.luma[,5:ncol(CNA.freqs.luma)], 1, function(x) (
    length(which(x==1))/ncol(CNA.freqs.luma[,5:ncol(CNA.freqs.luma)]))*100)
CNA.freqs.luma <- CNA.freqs.luma[c("gene","chr","start","end","freqloss.LumA","freqgain.LumA")]

# prep lumb files
ascat.df.lumb <- do.call(rbind, lapply(ascat.list.lumb, function(x) x$CNA))
ascat.df.lumb <- t(ascat.df.lumb)
ascat.df.lumb <- cbind(ascat.list.lumb[[1]][c("gene","chr","start","end")],ascat.df.lumb)
# get freqs
CNA.freqs.lumb <- ascat.df.lumb
CNA.freqs.lumb$freqloss.LumB <- apply(
  CNA.freqs.lumb[,5:ncol(CNA.freqs.lumb)], 1, function(x) (
    length(which(x==-1))/ncol(CNA.freqs.lumb[,5:ncol(CNA.freqs.lumb)]))*-100)
CNA.freqs.lumb$freqgain.LumB <- apply(
  CNA.freqs.lumb[,5:ncol(CNA.freqs.lumb)], 1, function(x) (
    length(which(x==1))/ncol(CNA.freqs.lumb[,5:ncol(CNA.freqs.lumb)]))*100)
CNA.freqs.lumb <- CNA.freqs.lumb[c("gene","chr","start","end","freqloss.LumB","freqgain.LumB")]

#######################################################################
# prepare and save files
#######################################################################

# exclude rows with NA
ascat.df.her2e <- ascat.df.her2e[complete.cases(ascat.df.her2e), ]
ascat.df.luma <- ascat.df.luma[complete.cases(ascat.df.luma), ]
ascat.df.lumb <- ascat.df.lumb[complete.cases(ascat.df.lumb), ]

# which genes are in common
common.genes <- intersect(intersect(ascat.df.her2e$gene, ascat.df.luma$gene), 
                          ascat.df.lumb$gene)

# only include common genes in ASCAT
ascat.df.her2e <- ascat.df.her2e[ascat.df.her2e$gene %in% common.genes, ]
ascat.df.luma <- ascat.df.luma[ascat.df.luma$gene %in% common.genes, ]
ascat.df.lumb <- ascat.df.lumb[ascat.df.lumb$gene %in% common.genes, ]

# only include common genes in frequencies
CNA.freqs.her2e <- CNA.freqs.her2e[CNA.freqs.her2e$gene %in% common.genes, ]
CNA.freqs.luma <- CNA.freqs.luma[CNA.freqs.luma$gene %in% common.genes, ]
CNA.freqs.lumb <- CNA.freqs.lumb[CNA.freqs.lumb$gene %in% common.genes, ]

# save
subtype.samples <- list("her2e"=colnames(ascat.df.her2e)[5:ncol(ascat.df.her2e)],
                        "LumA"=colnames(ascat.df.luma)[5:ncol(ascat.df.luma)],
                        "LumB"=colnames(ascat.df.lumb)[5:ncol(ascat.df.lumb)])
# make into one df
ascat.df.all <- merge(
  merge(ascat.df.her2e,ascat.df.luma[c(1,5:ncol(ascat.df.luma))],by="gene"),
  ascat.df.lumb[c(1,5:ncol(ascat.df.lumb))],by="gene")
ascat.save.file <- list("ascat.df.all" = ascat.df.all,
                        "subtype.samples" = subtype.samples)
save(ascat.save.file,file=outfile.1)

# make into one df
CNA.freqs.all <- merge(
  merge(CNA.freqs.her2e,CNA.freqs.luma[c("gene","freqloss.LumA","freqgain.LumA")],by="gene"),
  CNA.freqs.lumb[c("gene","freqloss.LumB","freqgain.LumB")], by="gene")
save(CNA.freqs.all,file=outfile.2)
