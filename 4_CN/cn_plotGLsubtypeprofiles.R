# Script: Plot gain/loss frequencies over the whole genome for all subtypes

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
gene.freqs.file <- "data/COMBINED/4_CN/processed/CNA_genefreqs.RData"
plot.file <- "output/plots/4_CN/COMBINED_HER2n_GLsubtypeprofiles.pdf"
#-------------------------------------------------------------------------------
# storing objects 
plot.list <- list() # object to store plots; note: if the output is not in string format use capture.output()
txt.out <- c() # object to store text output

###############################################################################
# loading data
###############################################################################

# gene level data
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

###############################################################################
# calc. gain and loss frequencies
###############################################################################

her2e.cna <- scanb.cna
luma.cna <- basis.cna[names(basis.cna) %in% c(c("gene","chr","centerPos","Genome_pos",luma.ids))]
lumb.cna <- basis.cna[c("gene","chr","centerPos","Genome_pos",lumb.ids)]

# check that genes are in same order in all dfs
# gene genome_pos freq.gain.her2e freq.loss.her2e freq.gain.luma freq.loss.luma freq.gain.lumb freq.loss.lumb
cna.list <- list("her2e"=her2e.cna,
     "luma"=luma.cna,
     "lumb"=lumb.cna)

for (i in names(cna.list)) {
  cna.df <- list("her2e"=her2e.cna,
                 "luma"=luma.cna,
                 "lumb"=lumb.cna)[[i]]
  freqs <- apply(cna.df, 1, function(y) {
    dat <- as.numeric(as.vector(y)[5:length(y)])
    freq.gain <- (length(dat[dat>0])/length(dat))*100
    freq.loss <- (length(dat[dat<0])/length(dat))*100
    return(c("gain"=freq.gain,"loss"=freq.loss))
  })
  freqs <- as.data.frame(freqs)
  names(freqs) <- cna.df$gene
  rownames(freqs) <- c(paste("freq.gain.",i,sep=""),
                       paste("freq.loss.",i,sep=""))
  #View(head(freqs))
  assign(paste(i,".freqs",sep=""), freqs)
}

# check that the gene/col order is the sample before rbdining
all(colnames(her2e.freqs)==colnames(luma.freqs)) # yes
all(colnames(her2e.freqs)==colnames(lumb.freqs)) # yes

# rbind
gene.freqs <- t(rbind(her2e.freqs,luma.freqs,lumb.freqs)) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene")
# add gene anno
gene.freqs$chr <- her2e.cna$chr[match(
  gene.freqs$gene,her2e.cna$gene)]
gene.freqs$centerPos <- her2e.cna$centerPos[match(
  gene.freqs$gene,her2e.cna$gene)]
gene.freqs$Genome_pos <- her2e.cna$Genome_pos[match(
  gene.freqs$gene,her2e.cna$gene)]

#gene.freqs <- gene.freqs[order(gene.freqs$Genome_pos),]
# convert loss to negative values
gene.freqs[c("freq.loss.her2e",
             "freq.loss.luma",
             "freq.loss.lumb")] <- lapply(gene.freqs[c("freq.loss.her2e",
                                                       "freq.loss.luma",
                                                       "freq.loss.lumb")], 
                                          FUN = function(x) {x*-1})

#View(head(gene.freqs))

# chr lengths file
chr.lengths <- as.data.frame(read.table(file = "data/BASIS/4_CN/raw/GRCh38_EBV.chrom.sizes.tsv", sep = '\t', header = FALSE))[1:23,] %>% 
    dplyr::rename(Chr=V1,length=V2) %>% 
  mutate(genome = cumsum(as.numeric(length))) %>% # dont shift by 1 here
  filter(Chr != "chrX")

#chr.lengths$Chr <- as.numeric(gsub("X",23,gsub('^.{3}','',chr.lengths$Chr)))
chr.lengths$Chr <- as.numeric(gsub('^.{3}','',chr.lengths$Chr))

# vector that has always the middle value for each chromosome (for plotting)
chr.lengths <- chr.lengths %>% mutate(chrbreaks = genome - length/2)

#View(cn.data)

save(gene.freqs,file=gene.freqs.file)
###############################################################################
# plot: all profiles without points
# The genome position is wrong
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

###############################################################################
# signif gene data
###############################################################################

# prep dat
gene.test.df <- loadRData("data/COMBINED/4_CN/processed/CNA_genelevel.RData")
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


###############################################################################
# plot: all profiles with points
###############################################################################

#pdf(file = paste(output.path,cohort,"_GLsubtypeprofiles.pdf", sep =""), height = 21.0, width = 72.0)

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
  geom_point(aes(x = genes.AG$Genome_pos, y = genes.AG$y), size=12) +
  geom_point(aes(x = genes.AL$Genome_pos, y = genes.AL$y), size=12) +
  geom_point(aes(x = genes.BG$Genome_pos, y = genes.BG$y), size=12, colour="red") +
  geom_point(aes(x = genes.BL$Genome_pos, y = genes.BL$y), size=12, colour="red") +
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

###############################################################################
# plot: separate profiles with points
###############################################################################

# luma

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
    x = cn.data[which(!is.na(cn.data$freq.gain.her2e)),]$Genome_pos, 
    y = cn.data[!is.na(cn.data$freq.gain.her2e),]$freq.gain.her2e, 
    color = "HER2E"),size=4) + 
  geom_line(aes(
    x = cn.data[which(!is.na(cn.data$freq.loss.her2e)),]$Genome_pos, 
    y = cn.data[!is.na(cn.data$freq.loss.her2e),]$freq.loss.her2e, 
    color = "HER2E"),size=4) + 
  scale_colour_manual(name="Subtype", values = c("HER2E"="#d334eb", "LUMA"="#2176d5", "LUMB"="#34c6eb")) + 
  geom_point(aes(x = genes.AG$Genome_pos, y = genes.AG$y), size=12) +
  geom_point(aes(x = genes.AL$Genome_pos, y = genes.AL$y), size=12) +
  geom_vline(xintercept = chr.lengths$genome[-length(chr.lengths$genome)], linetype="dashed",size=1) +
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


# lumb
plot <- ggplot() +  
  ggtitle("Genome-wide frequency of gain/loss CN alterations") +
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
  geom_point(aes(x = genes.BG$Genome_pos, y = genes.BG$y), size=12) +
  geom_point(aes(x = genes.BL$Genome_pos, y = genes.BL$y), size=12) +
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
