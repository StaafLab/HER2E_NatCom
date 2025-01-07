# Script: Plotting genome wide CNA Frequencies in SCAN-B 
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
pacman::p_load(IRanges,
               GenomicFeatures,
               TxDb.Hsapiens.UCSC.hg38.knownGene,
               TxDb.Hsapiens.UCSC.hg19.knownGene,
               org.Hs.eg.db,
               Repitools)
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
infile.4 <- "./data/BASIS/4_CN/raw/GRCh38_EBV.chrom.sizes.tsv"
#infile.5 <- "./data/Parameters/color_palette.RData"
infile.6 <- "./data/SCANB/4_CN/processed/CNA_GLFreqs_all.RData"
infile.7 <- "./data/SCANB/4_CN/processed/CNA_genetest.RData"
# output paths
#outfile.1 <- ""
plot.file <- paste0(output.path,cohort,"_CNAprofile.pdf")
#txt.file <- paste0(output.path,cohort,"_i.txt")
#-------------------
# storing objects 
#plot.list <- list() # object to store plots
#plot.parameters <- list() # object to store parameters to plot base R plots again later
#txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

color.palette <- c(LumA = "#0f1bbd", LumB = "#09d3e6", her2e ="#b609e6")

# load her2e ids
#her2e.ids <- unname(unlist(loadRData(infile.1)["ERpHER2n_her2e"]))

# load
gl.freqs <- loadRData(infile.6)

#basis.anno <- loadRData(infile.3)
#basis.anno <- basis.anno[basis.anno$ClinicalGroup == "ERposHER2neg" & basis.anno$PAM50_AIMS %in% c("LumA","LumB"),c("sample_name","PAM50_AIMS")]

#######################################################################
# add genome center position for genes
#######################################################################

# Read data from file
chr.lengths <- as.data.frame(read.table(file = infile.4, sep = '\t', header = FALSE))[1:23, ]

# Rename columns
names(chr.lengths) <- c("Chr", "length")

# Create 'genome' column
chr.lengths$genome <- cumsum(as.numeric(chr.lengths$length)) # lag by 1 position (cause I have to add the length of the previous chr to the probe positions (0 for chr1 probes))

# Shift 'genome' column by one row
chr.lengths$genome <- c(0, chr.lengths$genome[-nrow(chr.lengths)])

chr.lengths[chr.lengths$Chr=="chrX",]$Chr <- "chr23"
chr.lengths$Chr <- as.numeric(gsub('^.{3}','',chr.lengths$Chr))

# add center pos
gl.freqs[c("chr","start","end")] <- lapply(gl.freqs[c("chr","start","end")], as.numeric)
gl.freqs$centerPos <- (gl.freqs$start + gl.freqs$end) / 2

# convert to genome position
gl.freqs$genome <- chr.lengths$genome[match(gl.freqs$chr,chr.lengths$Chr)]

# Add a column with the genome position of each probe
gl.freqs$Genome_pos <- gl.freqs$centerPos + gl.freqs$genome

###############################################################################
# signif gene data
###############################################################################

# exclucde chr 23
gl.freqs <- gl.freqs[gl.freqs$chr!=23,]

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

# Pick higher frequency between freq.gain.luma and freq.gain.her2e
genes.AG$y <- with(genes.AG,ifelse(freqgain.LumA > freqgain.her2e,
                                   freqgain.LumA, freqgain.her2e))
                   
genes.AL$y <- with(genes.AL,ifelse(freqloss.LumA < freqloss.her2e,
                                   freqloss.LumA, freqloss.her2e))

genes.BG$y <- with(genes.BG,ifelse(freqgain.LumB > freqgain.her2e,
                                   freqgain.LumB, freqgain.her2e))

genes.BL$y <- with(genes.BL,ifelse(freqloss.LumB < freqloss.her2e,
                                   freqloss.LumB, freqloss.her2e))

###############################################################################
# plot: all profiles with points 
###############################################################################

pdf(file = plot.file, height = 21.0, width = 72.0)

plot <- ggplot() +  
  ggtitle("Genome-wide frequency of gain/loss CN alterations") +
  geom_line(aes(
    x = gl.freqs$Genome_pos, 
    y = gl.freqs$freqgain.LumA, 
    color = "LumA"),size=4) + 
  geom_line(aes(
    x = gl.freqs$Genome_pos, 
    y = gl.freqs$freqloss.LumA, 
    color = "LumA"),size=4) + 
  geom_line(aes(
    x = gl.freqs$Genome_pos, 
    y = gl.freqs$freqgain.LumB, 
    color = "LumB"),size=4) + 
  geom_line(aes(
    x = gl.freqs$Genome_pos, 
    y = gl.freqs$freqloss.LumB,  
    color = "LumB"),size=4) + 
  geom_line(aes(
    x = gl.freqs$Genome_pos, 
    y = gl.freqs$freqgain.her2e,  
    color = "her2e"),size=4) + 
  geom_line(aes(
    x = gl.freqs$Genome_pos, 
    y = gl.freqs$freqloss.her2e,  
    color = "her2e"),size=4) + 
  scale_colour_manual(name="Subtype", values = color.palette) + 
  geom_point(aes(x = genes.AG$Genome_pos, y = genes.AG$y), size=12) +
  geom_point(aes(x = genes.AL$Genome_pos, y = genes.AL$y), size=12) +
  geom_point(aes(x = genes.BG$Genome_pos, y = genes.BG$y), size=12, colour="red") +
  geom_point(aes(x = genes.BL$Genome_pos, y = genes.BL$y), size=12, colour="red") +
  geom_vline(xintercept = chr.lengths$genome[1:22],
             linetype="dashed",size=1) + 
  scale_x_continuous(name="Genome position (chromosome)",
                     breaks=chr.lengths$genome[1:22], 
                     labels=as.character(1:22),
                     limits = c(0,max(chr.lengths$genome[1:22])), #+50000000
                     expand = c(0, 0)) +
  scale_y_continuous(name="Alteration frequency (%)",
                     breaks=c(seq(-100,100,25)),
                     labels=c(100,75,50,25,0,25,50,75,100),
                     expand = c(0, 0),
                     limits = c(-100,100)) +
  theme_bw() +
  theme(text=element_text(size=30),
        legend.title = element_blank(),
        axis.title.y = element_text(vjust = 0.5),
        legend.position = c(0.97, 0.95),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",linewidth=2),
        axis.ticks = element_line(colour = "black", linewidth = 2),
        axis.ticks.length=unit(0.5, "cm")) + #legend.position = "none") +
  annotate(x=min(gl.freqs$Genome_pos)+30000000,
           y=c(-50,50), label=c("Loss","Gain"), 
           geom="text", angle=90, hjust=0.5, 
           size=9, colour=c("black","black")) 
print(plot)
#dev.off()

###############################################################################
# plot: all profiles with points specific to her2e
###############################################################################

# subset signif genes
genes.all.gain <- gl.freqs[gl.freqs$gene %in% 
                             intersect(genes.AG$gene,genes.BG$gene), ]
genes.all.loss <- gl.freqs[gl.freqs$gene %in% 
                             intersect(genes.AL$gene,genes.BL$gene), ]

# Pick higher frequency between freq.gain.luma and freq.gain.her2e
genes.all.gain$y <- apply(genes.all.gain[, c("freqgain.her2e", "freqgain.LumA", "freqgain.LumB")], 1, max)
genes.all.loss$y <- apply(genes.all.loss[, c("freqloss.her2e", "freqloss.LumA", "freqloss.LumB")], 1, min)

plot <- ggplot() +  
  ggtitle("Genome-wide frequency of gain/loss CN alterations") +
  geom_line(aes(
    x = gl.freqs$Genome_pos, 
    y = gl.freqs$freqgain.LumA, 
    color = "LumA"),size=4) + 
  geom_line(aes(
    x = gl.freqs$Genome_pos, 
    y = gl.freqs$freqloss.LumA, 
    color = "LumA"),size=4) + 
  geom_line(aes(
    x = gl.freqs$Genome_pos, 
    y = gl.freqs$freqgain.LumB, 
    color = "LumB"),size=4) + 
  geom_line(aes(
    x = gl.freqs$Genome_pos, 
    y = gl.freqs$freqloss.LumB,  
    color = "LumB"),size=4) + 
  geom_line(aes(
    x = gl.freqs$Genome_pos, 
    y = gl.freqs$freqgain.her2e,  
    color = "her2e"),size=4) + 
  geom_line(aes(
    x = gl.freqs$Genome_pos, 
    y = gl.freqs$freqloss.her2e,  
    color = "her2e"),size=4) + 
  scale_colour_manual(name="Subtype", values = color.palette) + 
  geom_point(aes(x = genes.all.gain$Genome_pos, y = genes.all.gain$y), size=12) +
  geom_point(aes(x = genes.all.loss$Genome_pos, y = genes.all.loss$y), size=12) +
  geom_vline(xintercept = chr.lengths$genome[-length(chr.lengths$genome)],
             linetype="dashed",size=1) + 
  scale_x_continuous(name="Genome position (chromosome)",
                     breaks=chr.lengths$genome, 
                     labels=as.character(1:23),
                     limits = c(0,max(chr.lengths$genome)), #+50000000
                     expand = c(0, 0)) +
  scale_y_continuous(name="Alteration frequency (%)",
                     breaks=c(seq(-100,100,25)),
                     labels=c(100,75,50,25,0,25,50,75,100),
                     expand = c(0, 0),
                     limits = c(-100,100)) +
  theme_bw() +
  theme(text=element_text(size=30),
        legend.title = element_blank(),
        axis.title.y = element_text(vjust = 0.5),
        legend.position = c(0.97, 0.95),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",linewidth=2),
        axis.ticks = element_line(colour = "black", linewidth = 2),
        axis.ticks.length=unit(0.5, "cm")) + #legend.position = "none") +
  annotate(x=min(gl.freqs$Genome_pos)+30000000,
           y=c(-50,50), label=c("Loss","Gain"), 
           geom="text", angle=90, hjust=0.5, 
           size=9, colour=c("black","black")) 
print(plot)
dev.off()