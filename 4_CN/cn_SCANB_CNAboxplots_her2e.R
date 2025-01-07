# Script: Plotting CNA boxplots in SCAN-B/BASIS
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
#infile.1 <- "./data/SCANB/0_GroupSamples/ERpHER2n_sampleIDs.RData"
#infile.2 <- "./data/SCANB/4_CN/processed/CNA_GLFreqs_all.RData"
#infile.3 <- "./data/BASIS/1_clinical/raw/Summarized_Annotations_BASIS.RData"
#infile.5 <- "./data/Parameters/color_palette.RData"
infile.6 <- "./data/SCANB/4_CN/processed/CNA_GLFreqs_all.RData"
infile.7 <- "./data/SCANB/4_CN/processed/CNA_genetest.RData"
infile.8 <- "./data/SCANB/4_CN/processed/CNA_genelevel_all.RData"
# output paths
#outfile.1 <- ""
plot.file <- paste0(output.path,cohort,"_CNAboxplots.pdf")
txt.file <- paste0(output.path,cohort,"_CNAboxplots.txt")
#-------------------
# storing objects 
#plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

color.palette <- c(her2e ="#b609e6", LumA = "#0f1bbd", LumB = "#09d3e6")

# load her2e ids
#her2e.ids <- unname(unlist(loadRData(infile.1)["ERpHER2n_her2e"]))

# load
gl.freqs <- loadRData(infile.6)

# exclucde chr 23
gl.freqs <- gl.freqs[gl.freqs$chr!=23,]

###############################################################################
# signif gene data
###############################################################################

# prep dat
gene.test.df <- loadRData(infile.7)
gene.test.df <- gene.test.df[gene.test.df$gene %in% gl.freqs$gene,]

# subset signif genes
genes.AG <- gl.freqs[gl.freqs$gene %in% 
                       gene.test.df[gene.test.df$LumA.Gain.padj <= 0.05, ]$gene, ]

genes.AL <- gl.freqs[gl.freqs$gene %in% 
                       gene.test.df[gene.test.df$LumA.Loss.padj <= 0.05, ]$gene, ]

genes.BG <- gl.freqs[gl.freqs$gene %in% 
                       gene.test.df[gene.test.df$LumB.Gain.padj <= 0.05, ]$gene, ]

genes.BL <- gl.freqs[gl.freqs$gene %in% 
                       gene.test.df[gene.test.df$LumB.Loss.padj <= 0.05, ]$gene, ]

###############################################################################
# plot: all profiles with points specific to her2e
###############################################################################

# # subset signif genes
genes.all.gain <- gl.freqs[gl.freqs$gene %in%
                             intersect(genes.AG$gene,genes.BG$gene), ]
genes.all.loss <- gl.freqs[gl.freqs$gene %in%
                             intersect(genes.AL$gene,genes.BL$gene), ]
# 
# # Pick higher frequency between freq.gain.luma and freq.gain.her2e
# genes.all.gain$y <- apply(genes.all.gain[, c("freqgain.her2e", "freqgain.LumA", "freqgain.LumB")], 1, max)
# genes.all.loss$y <- apply(genes.all.loss[, c("freqloss.her2e", "freqloss.LumA", "freqloss.LumB")], 1, min)
# 
# ####### ADAPT 

################################################################################
# CN boxplot freq differences of significant genes
################################################################################

# 2 data objects with luma and lumb
genes.AG.diff <- abs(genes.AG$freqgain.her2e - genes.AG$freqgain.LumA)
genes.AL.diff <- abs(genes.AL$freqloss.her2e - genes.AL$freqloss.LumA)
luma.dat <- c(genes.AG.diff,genes.AL.diff)
genes.BG.diff <- abs(genes.BG$freqgain.her2e - genes.BG$freqgain.LumB)
genes.BL.diff <- abs(genes.BL$freqloss.her2e - genes.BL$freqloss.LumB)
lumb.dat <- c(genes.BG.diff,genes.BL.diff)

## save source dat
# Create the LumA and LumB dataframes as before
# LumA genes - Gain and Loss differences
# lumA_gain_df <- data.frame(Gene = genes.AG$gene, LumA.diff = genes.AG.diff)
# lumA_loss_df <- data.frame(Gene = genes.AL$gene, LumA.diff = genes.AL.diff)
# # LumB genes - Gain and Loss differences
# lumB_gain_df <- data.frame(Gene = genes.BG$gene, LumB.diff = genes.BG.diff)
# lumB_loss_df <- data.frame(Gene = genes.BL$gene, LumB.diff = genes.BL.diff)
# # Merge LumA gains and losses (preferring losses, if both exist)
# lumA_df <- merge(lumA_gain_df, lumA_loss_df, by = "Gene", all = TRUE)
# # If both gain and loss exist for LumA, prefer loss, else use whichever is available
# lumA_df$LumA.diff <- ifelse(!is.na(lumA_df$LumA.diff.y), lumA_df$LumA.diff.y, lumA_df$LumA.diff.x)
# # Select only the Gene and LumA.diff columns
# lumA_df <- lumA_df[, c("Gene", "LumA.diff")]
# # Add Group as LumA
# lumA_df$Group <- "LumA"
# # Merge LumB gains and losses (preferring losses, if both exist)
# lumB_df <- merge(lumB_gain_df, lumB_loss_df, by = "Gene", all = TRUE)
# # If both gain and loss exist for LumB, prefer loss, else use whichever is available
# lumB_df$LumB.diff <- ifelse(!is.na(lumB_df$LumB.diff.y), lumB_df$LumB.diff.y, lumB_df$LumB.diff.x)
# # Select only the Gene and LumB.diff columns
# lumB_df <- lumB_df[, c("Gene", "LumB.diff")]
# # Add Group as LumB
# lumB_df$Group <- "LumB"
# # Merge LumA and LumB dataframes by Gene, keeping all genes (use all = TRUE)
# final_df <- merge(lumA_df, lumB_df, by = "Gene", all = TRUE)
# # Create a new combined Group column
# final_df$Group <- ifelse(!is.na(final_df$LumA.diff) & !is.na(final_df$LumB.diff), "both", 
#                          ifelse(!is.na(final_df$LumA.diff), "LumA", "LumB"))
# # Select only relevant columns
# final_df <- final_df[, c("Gene", "LumA.diff", "LumB.diff", "Group")]
# # View the final dataframe
# final_df
# save(final_df,file="./output/source_data/R_objects/Figure_4_genediff.RData")
# mean(final_df$LumA.diff,na.rm=TRUE)
# mean(final_df$LumB.diff,na.rm=TRUE)

##
mean(luma.dat)
mean(lumb.dat)

txt.out <- append(txt.out, c("\nGenes with significant difference in CNA Frequency compared to her2e-like\n",
                             "\n###########################################\n"))

txt.out <- append(txt.out, c("LumA Gain: n=",
                             paste0(nrow(genes.AG)," (",nrow(genes.AG)/nrow(gl.freqs)*100,"%)"),"\n",
                             "LumA Loss: n=",
                             paste0(nrow(genes.AL)," (",nrow(genes.AL)/nrow(gl.freqs)*100,"%)"),"\n",
                             "LumB Gain: n=",
                             paste0(nrow(genes.BG)," (",nrow(genes.BG)/nrow(gl.freqs)*100,"%)"),"\n",
                             "LumB Loss: n=",
                             paste0(nrow(genes.BL)," (",nrow(genes.BL)/nrow(gl.freqs)*100,"%)"),"\n",
                             "Both intersected Gain: n=",
                             paste0(nrow(genes.all.gain), " (",nrow(genes.all.gain)/nrow(gl.freqs)*100,"%)"),"\n",
                             "Both intersected Loss: n=",
                             paste0(nrow(genes.all.loss), " (",nrow(genes.all.loss)/nrow(gl.freqs)*100,"%)"),"\n",
                             "\n###########################################\n"))


plot.par <- list(
  data = list(LumA=luma.dat,LumB=lumb.dat), 
  col = color.palette[2:3], 
  names = names(color.palette[2:3]),
  ylab = "CNAFreq diff to her2e",
  main = "Abs. diff in CNA Freq to her2e")

plot.parameters <- append(plot.parameters, list(plot.par))

################################################################################
# CN boxplot of % genome altered
################################################################################

cna.df <- loadRData(infile.8)[[1]]
row.names(cna.df) <- cna.df$gene
sample.ids <- loadRData(infile.8)[[2]]

# sample # %altered
her2e.dat <- as.vector(apply(cna.df[sample.ids$her2e],2,function(x) {
  (sum(x!=0)/length(x))*100}))
luma.dat <-  as.vector(apply(cna.df[sample.ids$LumA],2,function(x) {
  (sum(x!=0)/length(x))*100}))
lumb.dat <-  as.vector(apply(cna.df[sample.ids$LumB],2,function(x) {
  (sum(x!=0)/length(x))*100}))

### source file export
# Function to calculate % genome altered for a given subtype
calc_percentage_altered <- function(samples, pam50_subtype) {
  altered_percent <- as.vector(apply(cna.df[samples], 2, function(x) {
    (sum(x != 0) / length(x)) * 100
  }))
  # Create a data frame for this subtype
  data.frame(
    SampleID = samples,
    PAM50 = pam50_subtype,
    GenomeAltered = altered_percent
  )
}
# Apply this function to each subtype and bind the results into a single data frame
her2e.dat.sf <- calc_percentage_altered(sample.ids$her2e, "HER2E")
luma.dat.sf <- calc_percentage_altered(sample.ids$LumA, "LumA")
lumb.dat.sf <- calc_percentage_altered(sample.ids$LumB, "LumB")
# Combine all the subtypes into a single data frame
genome_altered_df <- rbind(her2e.dat.sf, luma.dat.sf, lumb.dat.sf)
# View the final result
save(genome_altered_df, file= "./output/source_data/R_objects/Figure_4_genomealtered.RData")
#aggregate(GenomeAltered ~ PAM50, data = genome_altered_df, mean)
###

##
txt.out <- append(txt.out, c("\nGenome altered (%) compared to her2e\n",
                             "\n###########################################\n"))

txt.out <- append(txt.out, c("her2e: mean=", round(mean(her2e.dat),2),"\n",
                             "LumA: mean=", round(mean(luma.dat),2),"\n",
                             "LumB: mean=", round(mean(lumb.dat),2),"\n",
                             "\n###########################################\n"))

# mann whitney u tests
luma.res <- wilcox.test(her2e.dat, luma.dat) 
lumb.res <- wilcox.test(her2e.dat, lumb.dat)

txt.out <- append(txt.out, c(capture.output(luma.res), "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(lumb.res), "\n###########################################\n"))


plot.par <- list(
  data = list(her2e=her2e.dat,LumA=luma.dat,LumB=lumb.dat), 
  col = color.palette, 
  names = names(color.palette),
  ylab = "Genome altered (%)",
  main = "Genome altered")

plot.parameters <- append(plot.parameters, list(plot.par))

################################################################################
################################################################################

# save plots
pdf(file = plot.file, onefile = TRUE) 
par(mfrow = c(2, 2))
for (i in 1:length(plot.parameters)) {
  bp <- boxplot(plot.parameters[[i]]$data,
                col = plot.parameters[[i]]$col,
                names = plot.parameters[[i]]$names,
                ylab = plot.parameters[[i]]$ylab,
                main = plot.parameters[[i]]$main,
                ylim = c(0,100))
  axis(3,at=1:length(bp$n),labels=bp$n)
}
par(mfrow = c(1, 1))
dev.off()

# save text
writeLines(txt.out, txt.file)
