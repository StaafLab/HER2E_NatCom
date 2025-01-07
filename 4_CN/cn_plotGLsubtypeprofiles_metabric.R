# Script: METABRIC: Plot gain/loss frequencies over the whole genome for all subtypes

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

#-------------------------------------------------------------------------------
# set/create output directories
# for plots
output.path <- "output/plots/4_CN/"
dir.create(output.path)
# for data
data.path <- "data/METABRIC/4_CN/processed/"
dir.create(data.path)
#-------------------------------------------------------------------------------
# input & output file paths
# input
gene.cna <- "data/METABRIC/4_CN/processed/CNA_genelevel.RData"
gene.freqs.file <- "data/METABRIC/4_CN/processed/CN_gainloss_frequencies.RData"
# output
plot.file <- "output/plots/4_CN/METABRIC_HER2n_GLsubtypeprofiles.pdf"
#-------------------------------------------------------------------------------
# storing objects 
plot.list <- list() # object to store plots; note: if the output is not in string format use capture.output()
txt.out <- c() # object to store text output

###############################################################################
# loading data
###############################################################################

gene.freqs <- loadRData(gene.freqs.file)

#-----------------------------------------------------------------------------#

# signif genes
# prep dat
gene.test.df <- loadRData("data/METABRIC/4_CN/processed/Metabric_gene_cna_signif.RData")
genes.AG <- gene.freqs %>% 
  filter(gene %in% gene.test.df[gene.test.df$LumA.Gain.padj<=0.05,]$gene) %>% 
  mutate(y = ifelse(
    freq.gain.luma>freq.gain.her2e,freq.gain.luma,freq.gain.her2e)) #pick what is higher luma or her2e
genes.BG <- gene.freqs %>% 
  filter(gene %in% gene.test.df[gene.test.df$LumB.Gain.padj<=0.05,]$gene) %>% 
  mutate(y = ifelse(
    freq.gain.lumb>freq.gain.her2e,freq.gain.lumb,freq.gain.her2e))

genes.AL <- gene.freqs %>% 
  filter(gene %in% gene.test.df[gene.test.df$LumA.Loss.padj<=0.05,]$gene) %>% 
  mutate(y = ifelse(
    freq.loss.luma<freq.loss.her2e,freq.loss.luma,freq.loss.her2e))
genes.BL <- gene.freqs %>% 
  filter(gene %in% gene.test.df[gene.test.df$LumB.Loss.padj<=0.05,]$gene) %>% 
  mutate(y = ifelse(
    freq.loss.lumb<freq.loss.her2e,freq.loss.lumb,freq.loss.her2e))

#-----------------------------------------------------------------------------#


# add annotation: Biomart & Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",host="https://grch37.ensembl.org") # GRCh37

normal.chroms <- c(1:22)
my.symbols <- gene.freqs$gene

my.regions <- getBM(
  c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
  filters = c("hgnc_symbol", "chromosome_name"),
  values = list(hgnc_symbol=my.symbols, chromosome_name=normal.chroms),
  mart = ensembl)

# get data into the right format
gene.anno <- as.data.frame(my.regions) %>% 
  dplyr::rename(chr=chromosome_name) %>%
  mutate(pos=end_position-((end_position-start_position)/2)) %>% 
  dplyr::select(-c("start_position","end_position"))

# get chromosome lengths (take the last pos on each)
chr.lengths <- gene.anno %>% 
  group_by(chr) %>% 
  summarise(length = max(pos)) %>% 
  as.data.frame()
chr.lengths <- chr.lengths[with(chr.lengths, order(chr)),] %>% 
  mutate(genome = cumsum(length)) 
  #%>% mutate(chrbreaks = genome - length/2)

# add genome pos 
gene.freqs <- merge(dplyr::rename(gene.anno,"gene"="hgnc_symbol"),gene.freqs, by="gene") %>% 
  as.data.frame() %>% 
  rowwise() %>% 
  mutate(Genome_pos = ifelse(chr==1,pos,pos + sum(chr.lengths$length[1:chr-1]))) %>%
  relocate(Genome_pos, .after=pos)

# also add to signif gene data
genes.AG$Genome_pos <- gene.freqs$Genome_pos[match(
  genes.AG$gene,gene.freqs$gene)]
genes.AL$Genome_pos <- gene.freqs$Genome_pos[match(
  genes.AL$gene,gene.freqs$gene)]
genes.BG$Genome_pos <- gene.freqs$Genome_pos[match(
  genes.BG$gene,gene.freqs$gene)]
genes.BL$Genome_pos <- gene.freqs$Genome_pos[match(
  genes.BL$gene,gene.freqs$gene)]

###############################################################################
# plot: all profiles with points
###############################################################################
cn.data <- gene.freqs
pdf(file = plot.file, height = 21.0, width = 72.0)


plot <- ggplot() +  
  ggtitle("Genome-wide frequency of gain/loss CN alterations") +
  geom_line(aes(
    x = cn.data[which(!is.na(cn.data$freq.gain.luma)),]$Genome_pos, 
    y = cn.data[!is.na(cn.data$freq.gain.luma),]$freq.gain.luma, 
    color = "LUMA"),size=4) + 
  geom_line(aes(
    x = cn.data[which(!is.na(cn.data$freq.loss.luma)),]$Genome_pos, 
    y = cn.data[!is.na(cn.data$freq.loss.luma),]$freq.loss.luma, 
    color = "LUMA"),size=4) + 
  geom_line(aes(
    x = cn.data[which(!is.na(cn.data$freq.gain.lumb)),]$Genome_pos, 
    y = cn.data[!is.na(cn.data$freq.gain.lumb),]$freq.gain.lumb, 
    color = "LUMB"),size=4) + 
  geom_line(aes(
    x = cn.data[which(!is.na(cn.data$freq.loss.lumb)),]$Genome_pos, 
    y = cn.data[!is.na(cn.data$freq.loss.lumb),]$freq.loss.lumb, 
    color = "LUMB"),size=4) + 
  geom_line(aes(
    x = cn.data[which(!is.na(cn.data$freq.gain.her2e)),]$Genome_pos, 
    y = cn.data[!is.na(cn.data$freq.gain.her2e),]$freq.gain.her2e, 
    color = "HER2E"),size=4) + 
  geom_line(aes(
    x = cn.data[which(!is.na(cn.data$freq.loss.her2e)),]$Genome_pos, 
    y = cn.data[!is.na(cn.data$freq.loss.her2e),]$freq.loss.her2e, 
    color = "HER2E"),size=4) + 
  scale_colour_manual(name="Subtype", values = c("HER2E"="#d334eb", "LUMA"="#2176d5", "LUMB"="#34c6eb")) + 
  geom_point(aes(x = genes.AG$Genome_pos, y = genes.AG$y), size=9) +
  geom_point(aes(x = genes.AL$Genome_pos, y = genes.AL$y), size=9) +
  geom_point(aes(x = genes.BG$Genome_pos, y = genes.BG$y), size=9, colour="red") +
  geom_point(aes(x = genes.BL$Genome_pos, y = genes.BL$y), size=9, colour="red") +
  geom_vline(xintercept = chr.lengths$genome[-length(chr.lengths$genome)],
             linetype="dashed",size=1) +
  scale_x_continuous(name="Genome position (chromosome)",
                     breaks=chr.lengths$chrbreaks, 
                     labels=as.character(1:22),
                     limits = c(0,max(chr.lengths$genome)), #+50000000
                     expand = c(0, 0)) +
  scale_y_continuous(name="Alteration frequency (%)", # \n Loss          Gain", # \n Loss & G
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
  annotate(x=min(cn.data$Genome_pos)+30000000,
           y=c(-50,50), label=c("Loss","Gain"), 
           geom="text", angle=90, hjust=0.5, 
           size=9, colour=c("black","black")) 
print(plot)

dev.off()