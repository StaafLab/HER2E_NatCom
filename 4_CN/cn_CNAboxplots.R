# Script: plotting CN boxplot % of probes gain/loss; ploidy vs. PAM50; freq diff gain loss vs. pam50

#TODO: add stat testing?

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "COMBINED" #  

# set/create output directory for plots
output.path <- "output/plots/4_CN/"
dir.create(output.path)

# set/create output directory for data
data.path <- "data/SCANB/4_CN/processed/"
dir.create(data.path)

# plot
plot.list <- list() # object to store plots; note: if the output is not in string format use capture.output()
plot.file <- paste(output.path,cohort,"_HER2n_cnplots.pdf",sep = "")
txt.file <- paste(output.path,cohort,"_HER2n_cnplots_tests.txt", sep="")
txt.out <- c() # object to store text output

#packages
source("scripts/4_CN/src/cn_functions.R")
library(ggplot2)
library(tidyverse)
library(reshape2)
#library(data.table)
library(purrr)
library(readxl)

################################################################################
################################################################################

# 1. ploidy

# scanb
scanb.ploidy <- loadRData("./data/SCANB/4_CN/processed/HER2Esample_ploidy.RData") %>% 
  mutate(PAM50 = "Her2e") %>% dplyr::rename(Ploidy=ploidy)

# basis
basis.ploidy <- as.data.frame(read_excel("data/BASIS/4_CN/raw/Supplementary Table 5.Ploidy.AberrantCellFraction.090402015.v1.xlsx"))[1:2]
basis.anno <- loadRData("data/BASIS/1_clinical/raw/Summarized_Annotations_BASIS.RData") %>% 
  filter(ClinicalGroup == "ERposHER2neg" & PAM50_AIMS %in% c("LumA","LumB")) %>% 
  dplyr::select(c("sample_name","PAM50_AIMS")) %>% 
  dplyr::rename(Sample=sample_name,PAM50=PAM50_AIMS)
basis.ploidy <- merge(basis.anno,basis.ploidy,by="Sample")

all.ploidy <- rbind(scanb.ploidy, basis.ploidy)

p <- ggplot(all.ploidy, aes(x=PAM50,y=as.numeric(Ploidy),fill=PAM50)) +
  geom_boxplot(size=2.5, outlier.size = 7) +
  ylab("Tumor ploidy") +
  xlab("PAM50 subtype") +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",linewidth=2),
        axis.ticks = element_line(colour = "black", linewidth = 2),
        axis.ticks.length=unit(0.5, "cm")) +
  scale_fill_manual(values=setNames(c("#d334eb","#2176d5","#34c6eb"),
                                    c("Her2e","LumA","LumB"))) +  
  scale_y_continuous(breaks = scales::breaks_pretty(10))


print(p)
plot.list <- append(plot.list,list(p))#

################################################################################

# Statistics
# vars
her2e.data <- unname(unlist(subset(all.ploidy, PAM50 == "Her2e", select=c(Ploidy))))
for (subtype in c("LumA","LumB")) {
  # signature data for subtype
  comp.data <- unname(unlist(subset(all.ploidy, PAM50 == subtype, select=c(Ploidy))))
  
  test.type <- if(diff(range(unlist(her2e.data))) == 0 |
                  diff(range(unlist(comp.data))) == 0 ) {
    # all values identical for one or both -> cant be normally distrib.
    "non-parametric"
  } else if(shapiro.test(unlist(her2e.data))$p.value > 0.05 & 
            shapiro.test(unlist(comp.data))$p.value > 0.05)  {
    # parametric: both normally distributed
    "parametric"
  } else {     
    # at least one group not normally distributed
    "non-parametric"
  }
  
  if(test.type=="parametric") { 
    # equal variance check
    levenes.res <- var.test(unlist(her2e.data), unlist(comp.data), alternative = "two.sided")
    # t.test
    if (!is.nan(levenes.res$p.value) & levenes.res$p.value <= 0.05) {
      res <- t.test(her2e.data, comp.data, var.equal = FALSE)
    } else if (is.nan(levenes.res$p.value)) {
      res <- t.test(her2e.data, comp.data, var.equal = TRUE)
    } else { res <- t.test(her2e.data, comp.data, var.equal = TRUE) }
  } else { 
    # non-parametric
    # mann whitney u test
    res <- wilcox.test(unlist(her2e.data), unlist(comp.data), exact=TRUE)
  }
  # save result#
  #print(res)
  txt.out <- append(txt.out,c("Ploidy comparison",subtype,capture.output(res)))
}


################################################################################
# CN boxplot % of genes gain/loss 
################################################################################

# gene level data
scanb.cna <- loadRData(
  "data/SCANB/4_CN/processed/CNA_genelevel.RData")[["gainloss"]]
basis.cna <- loadRData("data/BASIS/4_CN/processed/CNA_genelevel_all.RData") %>% 
  dplyr::select(-c(chr,centerPos,Genome_pos))
basis.anno <- loadRData("data/BASIS/1_clinical/raw/Summarized_Annotations_BASIS.RData") %>% 
  filter(final.ER=="positive" & final.HER2=="negative") %>%
  dplyr::filter(PAM50_AIMS %in% c("LumA","LumB")) %>%
  dplyr::rename(PAM50=PAM50_AIMS) %>% 
  dplyr::select(sample_name,PAM50)

luma.ids <- basis.anno %>% filter(PAM50 == "LumA") %>% 
  filter(sample_name %in% colnames(basis.cna)) %>% 
  pull(sample_name)
lumb.ids <- basis.anno %>% filter(PAM50 == "LumB") %>% 
  filter(sample_name %in% colnames(basis.cna)) %>% 
  pull(sample_name)

# include only common genes
scanb.cna <- scanb.cna %>% drop_na() 
basis.cna <- basis.cna %>% drop_na()
common.genes <- intersect(scanb.cna$gene,basis.cna$gene)
scanb.cna <- scanb.cna %>% filter(gene %in% common.genes)
basis.cna <- basis.cna %>% filter(gene %in% common.genes)

# exclude not required columns
scanb.cna <- scanb.cna %>% 
  dplyr::select(-c("gene","chr","centerPos","Genome_pos"))
basis.cna <- basis.cna %>% 
  dplyr::select(-c("gene"))

# per sample: genes gain + genes loss / total genes = freq.altered
# final matrix: sample freq.altered PAM50

scanb.freq <- apply(scanb.cna, 2, function(x) {(length(x[x != 0])/length(x))*100}) %>% 
  as.data.frame() %>% 
  dplyr::rename("freq.altered"=".") %>% 
  rownames_to_column(var="sampleID") %>% 
  mutate(PAM50="HER2E")
basis.freq <- apply(basis.cna, 2, function(x) {(length(x[x != 0])/length(x))*100}) %>% 
  as.data.frame() %>% 
  dplyr::rename("freq.altered"=".") %>% 
  rownames_to_column(var="sampleID")
basis.freq$PAM50 <- basis.anno$PAM50[match(basis.freq$sampleID,basis.anno$sample_name)]

freqs <- rbind(basis.freq,scanb.freq)

# plot
p <- ggplot(freqs, aes(x=PAM50,y=as.numeric(freq.altered),fill=PAM50)) +
  geom_boxplot(size=2.5, outlier.size = 7) +
  ylab("Genome altered (%)") +
  xlab("PAM50 subtype") +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",linewidth=2),
        axis.ticks = element_line(colour = "black", linewidth = 2),
        axis.ticks.length=unit(0.5, "cm")) +
  scale_fill_manual(values=setNames(c("#d334eb","#2176d5","#34c6eb"),
                                    c("HER2E","LumA","LumB"))) +  
  scale_y_continuous(breaks = scales::breaks_pretty(10))


print(p)
plot.list <- append(plot.list,list(p))

################################################################################

# Statistics
# vars
her2e.data <- unname(unlist(subset(
  freqs, PAM50 == "HER2E", select=c(freq.altered))))
for (subtype in c("LumA","LumB")) {
  # signature data for subtype
  comp.data <- unname(unlist(subset(freqs, PAM50 == subtype, select=c(freq.altered))))
  
  test.type <- if(diff(range(unlist(her2e.data))) == 0 |
                  diff(range(unlist(comp.data))) == 0 ) {
    # all values identical for one or both -> cant be normally distrib.
    "non-parametric"
  } else if(shapiro.test(unlist(her2e.data))$p.value > 0.05 & 
            shapiro.test(unlist(comp.data))$p.value > 0.05)  {
    # parametric: both normally distributed
    "parametric"
  } else {     
    # at least one group not normally distributed
    "non-parametric"
  }
  
  if(test.type=="parametric") { 
    # equal variance check
    levenes.res <- var.test(unlist(her2e.data), unlist(comp.data), alternative = "two.sided")
    # t.test
    if (!is.nan(levenes.res$p.value) & levenes.res$p.value <= 0.05) {
      res <- t.test(her2e.data, comp.data, var.equal = FALSE)
    } else if (is.nan(levenes.res$p.value)) {
      res <- t.test(her2e.data, comp.data, var.equal = TRUE)
    } else { res <- t.test(her2e.data, comp.data, var.equal = TRUE) }
  } else { 
    # non-parametric
    # mann whitney u test
    res <- wilcox.test(unlist(her2e.data), unlist(comp.data), exact=TRUE)
  }
  # save result#
  #print(res)
  txt.out <- append(txt.out,c("Genome altered comparison ",subtype,capture.output(res)))
}

################################################################################
# CN boxplot freq differences of significant genes
################################################################################

# for the signfiicant genes, calc the difference in CNfreq and then visualize
genes <- loadRData("data/COMBINED/4_CN/processed/CNA_genelevel.RData")
#View(genes)
#colnames(genes)
luma.genes.loss <- genes %>% 
  filter(LumA.Loss.padj <= 0.05) %>% 
  pull(gene)
luma.genes.gain <- genes %>% 
  filter(LumA.Gain.padj <= 0.05) %>% 
  pull(gene)
lumb.genes.loss <- genes %>% 
  filter(LumB.Loss.padj <= 0.05) %>% 
  pull(gene)
lumb.genes.gain <- genes %>% 
  filter(LumB.Gain.padj <= 0.05) %>% 
  pull(gene)

# get the freq data and calc difference
gene.freqs <- loadRData("data/COMBINED/4_CN/processed/CNA_genefreqs.RData")

lumb.diff.loss <- gene.freqs %>% 
  filter(gene %in% lumb.genes.loss) %>% 
  mutate(diff = abs(freq.loss.her2e - freq.loss.lumb)) %>% 
  pull(diff)
lumb.diff.gain <- gene.freqs %>% 
  filter(gene %in% lumb.genes.gain) %>% 
  mutate(diff = abs(freq.gain.her2e - freq.gain.lumb)) %>% 
  pull(diff)
lumb.diff <- as.data.frame(c(lumb.diff.loss,lumb.diff.gain)) %>% 
  dplyr::rename("Freq.diff" = "c(lumb.diff.loss, lumb.diff.gain)") %>% 
  mutate(PAM50 = "LumB")

luma.diff.loss <- gene.freqs %>% 
  filter(gene %in% luma.genes.loss) %>% 
  mutate(diff = abs(freq.loss.her2e - freq.loss.luma)) %>% 
  pull(diff)
luma.diff.gain <- gene.freqs %>% 
  filter(gene %in% luma.genes.gain) %>% 
  mutate(diff = abs(freq.gain.her2e - freq.gain.luma)) %>% 
  pull(diff)
luma.diff <- as.data.frame(c(luma.diff.loss,luma.diff.gain)) %>% 
  dplyr::rename("Freq.diff" = "c(luma.diff.loss, luma.diff.gain)") %>% 
  mutate(PAM50 = "LumA")
freq.diffs <- rbind(lumb.diff,luma.diff)

# save the number of gain altered to text file
txt.out <- append(txt.out,
                  c("Number of genes signif. altered - Gain ",
                    paste("LumA: ",length(luma.genes.gain),"; LumB: ",length(lumb.genes.gain),sep="")))
txt.out <- append(txt.out,
                  c("Number of genes signif. altered - Loss ",
                    paste("LumA: ",length(luma.genes.loss),"; LumB: ",length(lumb.genes.loss),sep="")))


stats <- freq.diffs %>%
  group_by(PAM50) %>%
  summarise_each(funs(mean, median, sd), Freq.diff)
txt.out <- append(txt.out,
                  c("Signfic genes stats: Frequency diff in alterations",
                    capture.output(stats)))


# plot
p <- ggplot(freq.diffs, aes(x=PAM50,y=as.numeric(Freq.diff),fill=PAM50)) +
  geom_boxplot(size=2.5, outlier.size = 7) +
  ylab("Alteration frequency diff (%)") +
  xlab("PAM50 subtype") +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",linewidth=2),
        axis.ticks = element_line(colour = "black", linewidth = 2),
        axis.ticks.length = unit(0.5, "cm")) +
  scale_fill_manual(values = setNames(c("#2176d5","#34c6eb"),
                                    c("LumA","LumB"))) +  
  scale_y_continuous(breaks = scales::breaks_pretty(10))

#boxplot(Freq.diff~PAM50,data=freq.diffs,col=c("#2176d5","#34c6eb"))
#p <- recordPlot()
#plot.new()
plot.list <- append(plot.list,list(p))

########################################################################
########################################################################

# save plots
pdf(file = plot.file, onefile = TRUE)#, height = 5, width = 5)

for(i in 1:length(plot.list)) { 
   print(i)
   print(plot.list[[i]])
}

dev.off()

# save text output
writeLines(txt.out, txt.file)