# Script: HER2E project Mutational & Rearrangement signatures in SCAN-B
# Author: Lennart Hohmann
# Date: 19.02.2024
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
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readxl)
#-------------------
# set/create output directories
# for plots
output.path <- "./output/plots/3_genomic/"
dir.create(output.path)
# for data
data.path <- "./data/SCANB/3_genomic/processed/"
dir.create(data.path)
#-------------------
# input paths
infile.1 <- "./data/SCANB/3_genomic/raw/HER2_enriched_June23_ForJohan.xlsx"
infile.3 <- "./data/BASIS/1_clinical/raw/Summarized_Annotations_BASIS.RData"
infile.6 <- "./data/SCANB/3_genomic/raw/2024_02_11_snv_sig_refsig_errperc20_correct_exposures_combined_502.csv"
infile.9 <- "./data/SCANB/3_genomic/raw/2024_02_14_SCANB_ERpos_rearr_sig_502.csv"
# output paths
plot.file <- paste0(output.path,cohort,"_signatures.pdf")
txt.file <- paste0(output.path,cohort,"_signatures.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load new data
#######################################################################

qc.scanb <- as.data.frame(read_excel(infile.1, sheet = "Samples"))
pass.scanb.samples <- qc.scanb$Sample # all pass

# load palette
color.palette <- c(HER2E ="#b609e6", LumA = "#0f1bbd", LumB = "#09d3e6")

sign.mut <- read.table(infile.6, sep = ",", header = TRUE)
sign.mut <- sign.mut[c("Lund.tumour.id","SBS1","SBS2","SBS3","SBS5",
                     "SBS8","SBS13","SBS17","SBS18","SBS127", 
                     "unassigned","SBS6")]
names(sign.mut)[1] <- "Sample"
sign.mut$Sample <- gsub("\\..*", "", sign.mut$Sample)
rownames(sign.mut) <- sign.mut$Sample
sign.mut$Sample <- NULL

sign.rearr <- read.table(infile.9, sep = ",", header = TRUE)
sign.rearr <- sign.rearr[c("Lund.tumour.id","RefSigR2", "RefSigR4","RefSigR6a",
                         "RefSigR1","RefSigR5","RefSigR6b","RefSigR8","RefSigR3","RefSigR11",
                         "unassigned")]
names(sign.rearr)[1] <- "Sample"
sign.rearr$Sample <- gsub("\\..*", "", sign.rearr$Sample)
rownames(sign.rearr) <- sign.rearr$Sample
sign.rearr$Sample <- NULL

# add 6a and 6b
sign.rearr$RefSigR6 <- sign.rearr$RefSigR6a + sign.rearr$RefSigR6b
sign.rearr$RefSigR6a <- NULL
sign.rearr$RefSigR6b <- NULL

# which samples have <25 rearr snvs
sign.rearr.exclude <- sign.rearr
sign.rearr.exclude <- sign.rearr.exclude[pass.scanb.samples,]
sign.rearr.exclude$sum <-apply(sign.rearr.exclude,1,sum) #exclude and redo
sign.rearr.exclude <- row.names(sign.rearr.exclude[sign.rearr.exclude$sum < 25,])

# convert SCANB to proportion per sample
sign.rearr <- as.data.frame(apply(sign.rearr,1,function(x) (x/sum(x, na.rm=TRUE))))
sign.rearr[is.na(sign.rearr)] <- 0
sign.rearr <- sign.rearr[row.names(sign.rearr) != "unassigned",]

sign.mut <- as.data.frame(apply(sign.mut,1,function(x) (x/sum(x, na.rm=TRUE))))
sign.mut[is.na(sign.mut)] <- 0
sign.mut <- sign.mut[row.names(sign.mut) != "unassigned",]

# merge
sign.scanb <- as.data.frame(merge(t(sign.rearr), t(sign.mut), by = "row.names", all.x = TRUE))
names(sign.scanb)[1] <- "Sample"

# select HER2E samples
sign.scanb <- sign.scanb[sign.scanb$Sample %in% pass.scanb.samples,]

#######################################################################
# load and process BASIS data
#######################################################################

basis.anno <- loadRData(infile.3)
sign.basis <- basis.anno[basis.anno$ClinicalGroup == "ERposHER2neg" & 
                           basis.anno$PAM50_AIMS %in% c("LumA","LumB"),
                         c("PAM50_AIMS","RSproportions","MutSigProportions")]
# flatten
rs.df <- as.data.frame(sign.basis$RSproportions)
ms.df <- as.data.frame(sign.basis$MutSigProportions)
sign.basis <- cbind(sign.basis["PAM50_AIMS"], rs.df, ms.df)

# check which signatures are present in both cohorts
colnames(sign.scanb) <- gsub("RefSigR","RS", colnames(sign.scanb))
colnames(sign.scanb) <- gsub("SBS","S", colnames(sign.scanb))
common.sigs <- intersect(colnames(sign.scanb),colnames(sign.basis))

# exclude non-shared ones
rownames.her2e <- sign.scanb$Sample
sign.scanb <- sign.scanb[,common.sigs]
pam50.basis <- sign.basis$PAM50_AIMS
sign.basis <- sign.basis[,common.sigs]
sign.basis$PAM50 <- pam50.basis

# add subtype
sign.scanb$PAM50 <- "HER2E"

# put together in an object ready for plotting 
sign.dat <- rbind(sign.scanb,sign.basis)
#View(sign.dat)
row.names(sign.dat)[1:length(rownames.her2e)] <- rownames.her2e

# export to source file
sign.dat.sf <- sign.dat
sign.dat.sf$sampleID <- rownames(sign.dat.sf)
rownames(sign.dat.sf) <- NULL
sign.dat.sf <- sign.dat.sf[c("sampleID","PAM50","S1","S2","S5","S13","RS2")]
sign.dat.sf <- sign.dat.sf[!is.na(sign.dat.sf$PAM50),]
save(sign.dat.sf, file="./output/source_data/R_objects/Figure_4_sigs.RData")

#######################################################################
# plot
#######################################################################
sig="RS2"
for (sig in common.sigs) {
  
  txt.out <- append(txt.out, c("\n",sig,"\n",
                               "\n###########################################\n"))
  
  
  
  # sig data
  luma.dat <- sign.dat[sign.dat$PAM50=="LumA",sig]
  lumb.dat <- sign.dat[sign.dat$PAM50=="LumB",sig]
  # have to filter samples for rearr signatures below <25
  if (grepl("RS", sig)) {
    her2e.dat <- sign.dat[sign.dat$PAM50=="HER2E",]
    her2e.dat <- her2e.dat[!row.names(her2e.dat) %in% sign.rearr.exclude,sig]
  } else { her2e.dat <- sign.dat[sign.dat$PAM50=="HER2E",sig] }
  
  
  # statistics
  # summary statistics
  her2e.stats <- get_stats(her2e.dat)
  luma.stats <- get_stats(luma.dat)
  lumb.stats <- get_stats(lumb.dat)
  
  txt.out <- append(txt.out, c("HER2E\n",capture.output(her2e.stats), "\n",
                               "LumA\n",capture.output(luma.stats), "\n",
                               "LumB\n",capture.output(lumb.stats),
                               "\n###########################################\n"))
  
  # mann whitney u tests
  luma.res <- wilcox.test(her2e.dat, luma.dat)
  lumb.res <- wilcox.test(her2e.dat, lumb.dat)
  
  txt.out <- append(txt.out, c(capture.output(luma.res), "\n###########################################\n"))
  txt.out <- append(txt.out, c(capture.output(lumb.res), "\n###########################################\n"))
  
  # plot
  plot.par <- list(
    data = list(HER2E=her2e.dat,LumA=luma.dat,LumB=lumb.dat), 
    col = color.palette, 
    names = names(color.palette),
    ylab = "proportion",
    main = sig)
  plot.parameters <- append(plot.parameters, list(plot.par))

}

#######################################################################
#######################################################################

# save plots
pdf(file = plot.file, onefile = TRUE) 
par(mfrow = c(2, 2))
for (i in 1:length(plot.parameters)) {
  bp <- boxplot(plot.parameters[[i]]$data,
                col = plot.parameters[[i]]$col,
                names = plot.parameters[[i]]$names,
                ylab = plot.parameters[[i]]$ylab,
                main = plot.parameters[[i]]$main,
                ylim= c(0,1))
  axis(3,at=1:length(bp$n),labels=bp$n)
}
par(mfrow = c(1, 1))
dev.off()

# save text
writeLines(txt.out, txt.file)
