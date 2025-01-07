# Script: OoncotypeDX in ERp SCAN-B
# Author: Lennart Hohmann
# Date: 24.05.2024
#-------------------
# empty environment
rm(list=ls())
# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")
# cohort
cohort <- "SCANB"
#-------------------
# packages
source("scripts/4_CN/src/cn_functions.R")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(genefu,biomaRt,readxl)
#-------------------
# set/create output directories
# for plots
output.path <- "./output/plots/1_clinical/"
dir.create(output.path)
# for data
data.path <- "data/SCANB/1_clinical/processed/"
dir.create(data.path)
#-------------------
# input paths
infile.2 <- "./data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData"
infile.3 <- "./data/SCANB/2_transcriptomic/raw/genematrix_noNeg.Rdata"
infile.4 <- "../Project_Basal/data/SCANB/2_transcriptomic/processed/ERp_LogScaled_gex.RData"
infile.5 <- "../Project_Basal/data/SCANB/1_clinical/raw/SSP_SupplDataTable1.xlsx"
# output paths
plot.file <- paste0(output.path,cohort,"_oncoDX&ROR.pdf")
txt.file <- paste0(output.path,cohort,"_oncoDX&ROR.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots 
plot.parameters <- list() # object to store parameters to plot base R plots again later 
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# functions
#######################################################################

oldid_to_newid <- function(ids,
                            old_id="ensembl_gene_id",
                            new_id1="entrezgene_id",
                            new_id2="hgnc_symbol") { #entrezgene_id #hgnc_symbol
  # Connect to the Ensembl database
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  # Retrieve Hugo gene names using biomaRt
  res <- getBM(attributes = c(old_id, new_id1, new_id2),
               filters = old_id,
               values = ids,
               mart = mart,
               uniqueRows = TRUE)
  return(res)
}

get_stats <- function(vec) {
  mean <- mean(vec)
  median <- median(vec)
  sd <- sd(vec)
  return(c("mean"=mean, "median"=median, "sd"=sd))
}

#######################################################################
# load data
#######################################################################


# 
data(sig.oncotypedx)

# load annotation data and select 
anno <- loadRData(infile.2)
anno <- anno[anno$Follow.up.cohort == TRUE & anno$ER == "Positive" & anno$HER2 == "Negative", ] 
anno <- anno[anno$NCN.PAM50 %in% c("Her2","LumA","LumB"), ] 
# only include groups of interests (her2e, luma, lumb)

# lod ssp ror data
ror.dat <- read_excel(infile.5)
ror.dat <- ror.dat[ror.dat$Follow.up.cohort == TRUE & 
                     ror.dat$Sample %in% anno$Sample,
                   c("Sample","NCN.PAM50","NCN.ROR.risk.cat",
                     "NCN.ROR.binary.risk.cat","Size.mm",
                     "NCN.ROR.asT0","NCN.ROR.asT1")]
ror.dat <- ror.dat[!is.na(ror.dat$NCN.ROR.risk.cat), ]
ror.dat$NCN.ROR.comb <- ifelse(ror.dat$Size.mm <= 20,
                               ror.dat$NCN.ROR.asT0,
                               ror.dat$NCN.ROR.asT1)

# keep only samples present within both
anno <- anno[anno$Sample %in% ror.dat$Sample,]

# load and prep gene expr. data
gex <- as.data.frame(loadRData(infile.3))
gex <- as.data.frame(log2(gex + 1)) #log transform
#gex.scaled <- as.data.frame(loadRData(infile.4))
#gex.scaled[1:3,1:3]
#gex[1:3,(ncol(gex)-3):ncol(gex)]
# correct colnames
gex <- gex[anno$GEX.assay]
names(gex) <- anno$Sample[match(colnames(gex),anno$GEX.assay)]

#length(ror.dat$Sample)
#length(anno$Sample)
#ncol(gex)

#######################################################################
# calc risk scores
#######################################################################

# correct gene names
gex$Ensemble_ID <- gsub("\\..*$", "", rownames(gex))
rownames(gex) <- NULL
id.table <- oldid_to_newid(gex$Ensemble_ID, 
                           old_id="ensembl_gene_id")
gex$new.id <- id.table$entrezgene_id[match(gex$Ensemble_ID,id.table$ensembl_gene_id)]
gex <- gex[!duplicated(gex$new.id), ] # keep 1st row of each symbol
gex <- gex[!is.na(gex$new.id),]
rownames(gex) <- gex$new.id
gex$new.id <- NULL
gex$Ensemble_ID <- NULL

t.gex <- as.data.frame(t(gex))
t.gex <- t.gex[as.character(sig.oncotypedx$EntrezGene.ID)]
#View(t.gex)

# Matrix of annotations with at least one column named "EntrezGene.ID", dimnames being properly defined. 
# rename my genes
real.names <- colnames(t.gex)
colnames(t.gex) <- paste("G", seq(1, 21), sep="")
annot.dat <- data.frame(probe = colnames(t.gex), 
                        EntrezGene.ID = real.names)
rownames(annot.dat) <- colnames(t.gex)
res <- oncotypedx(data=t.gex, annot=annot.dat, 
                  do.mapping=TRUE,
                  do.scaling=TRUE)

#######################################################################
# summarize results: score 
#######################################################################
dat <- res$score

#### save source data ####
source_dat_path <- "./output/source_data/R_objects/"  
groups <- c("Her2", "LumA", "LumB")
results.OncoDX.score <- setNames(
  lapply(groups, function(group) {
    # Get the sample indices for the current group
    indices <- anno$Sample[anno$NCN.PAM50 == group]
    # Create the data vector and sample IDs
    list(data_vector = as.vector(dat[indices]),
         sample_ids = indices)
  }), 
  groups)
#saveRDS(results.OncoDX.score, paste0(source_dat_path,"OncoDX_score.RData")) # Save 
#### exported source data ####

# subtype data
her2.dat <- as.numeric(as.vector(
  dat[anno$Sample[anno$NCN.PAM50 =="Her2"]]))
luma.dat <- as.numeric(as.vector(
  dat[anno$Sample[anno$NCN.PAM50 =="LumA"]]))
lumb.dat <- as.numeric(as.vector(
  dat[anno$Sample[anno$NCN.PAM50 =="LumB"]]))

# summary statistics
her2.stats <- get_stats(her2.dat)
luma.stats <- get_stats(luma.dat)
lumb.stats <- get_stats(lumb.dat)

txt.out <- append(txt.out, c("Her2e\n",capture.output(her2.stats), "\n",
                             "LumA\n",capture.output(luma.stats), "\n",
                             "LumB\n",capture.output(lumb.stats),
                             "\n###########################################\n"))

# mann whitney u tests
luma.res <- wilcox.test(her2.dat, luma.dat)
lumb.res <- wilcox.test(her2.dat, lumb.dat)

txt.out <- append(txt.out, c(capture.output(luma.res), "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(lumb.res), "\n###########################################\n"))

# plot
plot.par <- list(
  data = list(LumA=luma.dat,LumB=lumb.dat,Her2e=her2.dat), 
  col = c(LumA = "#2176d5", LumB = "#34c6eb",HER2E = "#d334eb"), 
  names = names(c(LumA = "#2176d5", LumB = "#34c6eb",HER2E = "#d334eb")),
  ylab = "OncotypeDX scores")

plot.parameters <- append(plot.parameters, list(plot.par))
#plot.list <- append(plot.list, list(plot))


#######################################################################
# summarize results: risk
#######################################################################

dat <- res$risk

#### save source data ####
#source_dat_path <- "./output/source_data/R_objects/"  
#groups <- c("Her2", "LumA", "LumB")
results.OncoDX.group <- setNames(
  lapply(groups, function(group) {
    # Get the sample indices for the current group
    indices <- anno$Sample[anno$NCN.PAM50 == group]
    # Create the data vector and sample IDs
    list(data_vector = as.vector(dat[indices]),
         sample_ids = indices)
  }), 
  groups)
#saveRDS(results.OncoDX.group, paste0(source_dat_path,"OncoDX_group.RData")) # Save 
#### exported source data ####


# subtype data
# her2.dat <- as.numeric(as.vector(
#   dat[anno$Sample[anno$NCN.PAM50 =="Her2"]]))
# luma.dat <- as.numeric(as.vector(
#   dat[anno$Sample[anno$NCN.PAM50 =="LumA"]]))
# lumb.dat <- as.numeric(as.vector(
#   dat[anno$Sample[anno$NCN.PAM50 =="LumB"]]))

dat.an <- as.data.frame(dat)
rownames(anno) <- anno$Sample
dat.an <- merge(dat.an, anno[,c("NCN.PAM50"), drop = FALSE], by="row.names")

# sample ids for later comparison
dx.high <- dat.an[which(dat.an$dat==1),c("Row.names")]
dx.high.her2e <- dat.an[which(dat.an$dat==1 & dat.an$NCN.PAM50=="Her2"),c("Row.names")]

ct <- table(dat.an$NCN.PAM50,dat.an$dat)
print(ct)

res.luma <- fisher.test(ct[c("Her2","LumA"),])
res.lumb <- fisher.test(ct[c("Her2","LumB"),])

txt.out <- append(txt.out, 
                  c("\nOncotypeDX_riskgroup\n",
                    "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(res.luma), "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(res.lumb), "\n###########################################\n"))

variable <- c("low","intermediate","high")
her.counts <- unname(ct["Her2",])
luma.counts <- unname(ct["LumA",])
lumb.counts <- unname(ct["LumB",])
her.perc <- round((her.counts/sum(ct["Her2",]))*100,2)
luma.perc <- round((luma.counts/sum(ct["LumA",]))*100,2)
lumb.perc <- round((lumb.counts/sum(ct["LumB",]))*100,2)

ot.df <- as.data.frame(cbind(variable,
                              her.counts,her.perc,
                              luma.counts,luma.perc,
                              lumb.counts,lumb.perc))

#names(ot.df) <- x

#ot.df

txt.out <- append(txt.out, c(capture.output(ot.df), "\n###########################################\n"))


#######################################################################
# PROSIGNA ROR SCORE CONTINUE HERE
#######################################################################

#######################################################################
# summarize results: score 
#######################################################################
dat <- setNames(ror.dat$NCN.ROR.comb, ror.dat$Sample)

#### save source data ####
#source_dat_path <- "./output/source_data/R_objects/"  
#groups <- c("Her2", "LumA", "LumB")
results.ROR.score <- setNames(
  lapply(groups, function(group) {
    # Get the sample indices for the current group
    indices <- anno$Sample[anno$NCN.PAM50 == group]
    # Create the data vector and sample IDs
    list(data_vector = as.vector(dat[indices]),
         sample_ids = indices)
  }), 
  groups)
#saveRDS(results.ROR.score, paste0(source_dat_path,"ROR_score.RData")) # Save 
#### exported source data ####


# subtype data
her2.dat <- as.numeric(ror.dat[ror.dat$NCN.PAM50 =="Her2",]$NCN.ROR.comb)
luma.dat <- as.numeric(ror.dat[ror.dat$NCN.PAM50 =="LumA",]$NCN.ROR.comb)
lumb.dat <- as.numeric(ror.dat[ror.dat$NCN.PAM50 =="LumB",]$NCN.ROR.comb)

# summary statistics
her2.stats <- get_stats(her2.dat)
luma.stats <- get_stats(luma.dat)
lumb.stats <- get_stats(lumb.dat)

txt.out <- append(txt.out, 
                  c("\n###########################################\n",
                  "\n# ROR SCORE #\n",
                  "\n###########################################\n"))

txt.out <- append(txt.out, c("Her2e\n",capture.output(her2.stats), "\n",
                             "LumA\n",capture.output(luma.stats), "\n",
                             "LumB\n",capture.output(lumb.stats),
                             "\n###########################################\n"))

# mann whitney u tests
luma.res <- wilcox.test(her2.dat, luma.dat)
lumb.res <- wilcox.test(her2.dat, lumb.dat)

txt.out <- append(txt.out, c(capture.output(luma.res), "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(lumb.res), "\n###########################################\n"))

#######################################################################
# summarize results: risk
#######################################################################

dat <- setNames(ror.dat$NCN.ROR.risk.cat, ror.dat$Sample)

#### save source data ####
#source_dat_path <- "./output/source_data/R_objects/"  
#groups <- c("Her2", "LumA", "LumB")
results.ROR.group <- setNames(
  lapply(groups, function(group) {
    # Get the sample indices for the current group
    indices <- anno$Sample[anno$NCN.PAM50 == group]
    # Create the data vector and sample IDs
    list(data_vector = as.vector(dat[indices]),
         sample_ids = indices)
  }), 
  groups)
#saveRDS(results.ROR.group, paste0(source_dat_path,"ROR_group.RData")) # Save 
#### exported source data ####

# comp to oncodx class
her2e.ids <- anno$Sample[anno$NCN.PAM50=="Her2"]
ror.high <- ror.dat[which(ror.dat$NCN.ROR.risk.cat=="High"),][["Sample"]]
ror.high.her2e <- ror.dat[which(
  ror.dat$NCN.ROR.risk.cat=="High" & 
    ror.dat$NCN.PAM50=="Her2"),][["Sample"]]

# oncoDX high in both
common_high <- intersect(dx.high, ror.high)
percent_dxhigh <- (length(common_high) / length(dx.high)) * 100

#common_high <- intersect(ror.high, dx.high)
#percent_rorhigh <- (length(common_high) / length(ror.high)) * 100

# her2e high in both
num_common_ids <- sum(her2e.ids %in% common_high)
percent_her2e <- (num_common_ids / length(her2e.ids)) * 100

ct <- table(ror.dat$NCN.PAM50,ror.dat$NCN.ROR.risk.cat)

print(ct)

res.luma <- fisher.test(ct[c("Her2","LumA"),])
res.lumb <- fisher.test(ct[c("Her2","LumB"),])

txt.out <- append(txt.out, 
                  c("\nROR_riskgroup\n",
                    "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(res.luma), "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(res.lumb), "\n###########################################\n"))

variable <- c("high","intermediate","low")
her.counts <- unname(ct["Her2",])
luma.counts <- unname(ct["LumA",])
lumb.counts <- unname(ct["LumB",])
her.perc <- round((her.counts/sum(ct["Her2",]))*100,2)
luma.perc <- round((luma.counts/sum(ct["LumA",]))*100,2)
lumb.perc <- round((lumb.counts/sum(ct["LumB",]))*100,2)

ot.df <- as.data.frame(cbind(variable,
                             her.counts,her.perc,
                             luma.counts,luma.perc,
                             lumb.counts,lumb.perc))

#names(ot.df) <- x

#ot.df

txt.out <- append(txt.out, c(capture.output(ot.df), "\n###########################################\n"))


#######################################################################
#######################################################################

# save source data, make 1 df out of it
df1 <- data.frame(SampleID = unlist(lapply(results.OncoDX.group, function(x) x$sample_ids)),
                  OncoDX.group = unlist(lapply(results.OncoDX.group, function(x) x$data_vector)))

df2 <- data.frame(SampleID = unlist(lapply(results.OncoDX.score, function(x) x$sample_ids)),
                  OncoDX.score = unlist(lapply(results.OncoDX.score, function(x) x$data_vector)))

df3 <- data.frame(SampleID = unlist(lapply(results.ROR.group, function(x) x$sample_ids)),
                  ROR.group = unlist(lapply(results.ROR.group, function(x) x$data_vector)))

df4 <- data.frame(SampleID = unlist(lapply(results.ROR.score, function(x) x$sample_ids)),
                  ROR.score = unlist(lapply(results.ROR.score, function(x) x$data_vector)))

# Merge the data frames on SampleID
combined_df <- Reduce(function(x, y) merge(x, y, by = "SampleID", all = TRUE), 
                      list(df1, df2, df3, df4))

combined_df$PAM50 <- anno$NCN.PAM50[match(combined_df$SampleID,anno$Sample)] 
#table(combined_df$PAM50)
combined_df <- combined_df[, c("SampleID", "PAM50", 
                               "OncoDX.group", "OncoDX.score", 
                               "ROR.group", "ROR.score")]

#View(combined_df)
saveRDS(combined_df, paste0(source_dat_path,"Table_S3.RData")) # Save 


# save plots
pdf(file = plot.file, onefile = TRUE) 
#par(mfrow = c(2, 2))
for (i in 1:length(plot.parameters)) {
  bp <- boxplot(plot.parameters[[i]]$data,
                col = plot.parameters[[i]]$col,
                names = plot.parameters[[i]]$names,
                ylab = plot.parameters[[i]]$ylab)
  axis(3,at=1:length(bp$n),labels=bp$n)
}
par(mfrow = c(1, 1))
dev.off()

# save text
writeLines(txt.out, txt.file)
