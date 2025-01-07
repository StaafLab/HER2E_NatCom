# Script: Copy number alteration analyses in the Metabric cohort

#TODO: 

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "METABRIC"

# set/create output directory for plots
output.path <- "output/plots/4_CN/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/4_CN/processed/", sep="")
dir.create(data.path)
plot.file <- paste(output.path,cohort,"_HER2n_cnplots.pdf",sep = "")

txt.file <- paste(output.path,cohort,"_HER2n_cnplots_tests.txt", sep="")
txt.out <- c() # object to store text output

# store plots
plot.list <- list()

#packages
source("scripts/4_CN/src/cn_functions.R")
source("scripts/2_transcriptomic/src/tscr_functions.R")
library(ggplot2)
library(tidyverse)
library(matrixStats)
library(pheatmap)
library(VennDiagram)
library(readxl)
library(ggfortify)
library(janitor)
library(car)
library(ggsignif)
library(biomaRt)

#######################################################################
# 2. Cohort-specific data preprocessing including selection of  
# the clinical ER+Her2- subtyped samples
#######################################################################

# 2.1 load data
cn.data <- read.delim("./data/METABRIC/4_CN/raw/data_CNA.txt",header = TRUE, sep = "\t", dec = ".") %>% 
  dplyr::select(-c(Entrez_Gene_Id)) %>% 
  column_to_rownames(var = "Hugo_Symbol")

# annotation data
# load anno data
anno <- loadRData("data/METABRIC/1_clinical/processed/Merged_annotations_ERpHER2n.RData") %>% 
  dplyr::rename(sampleID=METABRIC_ID) %>% # rename to match SCANB variables
  mutate(sampleID = gsub("-",".",sampleID)) %>% 
  dplyr::select(sampleID,PAM50)
# anno <- loadRData("./data/METABRIC/1_clinical/raw/Merged_annotations.RData") %>%
#   as.data.frame(.) %>%
#   filter(PAM50 %in% c("LumA", "LumB", "Her2")) %>%
#   filter(grepl('ERpHER2n', ClinGroup)) %>%
#   dplyr::rename(sampleID=METABRIC_ID) %>%
#   dplyr::select(sampleID,PAM50) %>%
#   mutate(sampleID = gsub("-",".",sampleID)) # have to get identifiers in same format

cn.data <- cn.data  %>% 
  dplyr::select(any_of(anno$sampleID))

anno <- anno %>% 
  filter(sampleID %in% colnames(cn.data))

### quick FGFR4 check - put in manuscript
#prop.table(table(t(cn.data["FGFR4",c(anno[which(anno$PAM50=="Her2"),]$sampleID)])))*100
#prop.table(table(t(cn.data["FGFR4",c(anno[which(anno$PAM50=="LumA"),]$sampleID)])))*100
#prop.table(table(t(cn.data["FGFR4",c(anno[which(anno$PAM50=="LumB"),]$sampleID)])))*100
  
##########################################################################
# Part 1: general gain/loss, whole genome, freq. altered per group boxplot
##########################################################################
#head(cn.freq)

cn.freq <- apply(cn.data, 2, function(x) {(length(x[which(x != 0 & !is.na(x))])/length(x[!is.na(x)]))*100}) %>% 
  as.data.frame() %>% 
  dplyr::rename("freq.altered"=".") %>% 
  rownames_to_column(var="sampleID") 

# add to anno
cn.freq$PAM50 <- anno$PAM50[match(cn.freq$sampleID,anno$sampleID)]

# plot
p <- ggplot(cn.freq, aes(x=PAM50,y=as.numeric(freq.altered),fill=PAM50)) +
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
                                    c("Her2","LumA","LumB"))) +  
  scale_y_continuous(breaks = scales::breaks_pretty(10))


print(p)
plot.list <- append(plot.list,list(p))

################################################################################
freqs <- cn.freq
# Statistics
# vars
her2e.data <- unname(unlist(subset(
  freqs, PAM50 == "Her2", select=c(freq.altered))))
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
  print(res)
  txt.out <- append(txt.out,c("Genome altered comparison ",subtype,capture.output(res)))
}

################################################################################
# CN boxplot freq differences of significant genes
################################################################################
cn.data.test <- cn.data %>% na.omit()
#cn.data.test[1:5,1:5]
#row.names(cn.data.test)

her2e.data <- cn.data.test[,colnames(cn.data.test) %in% 
                            anno[anno$PAM50=="Her2",]$sampleID]
luma.data <- cn.data.test[,colnames(cn.data.test) %in% 
                            anno[anno$PAM50=="LumA",]$sampleID]
lumb.data <- cn.data.test[,colnames(cn.data.test) %in% 
                            anno[anno$PAM50=="LumB",]$sampleID]
#View(her2e.dat)

# save obj
# gene LumA.Gain.pval LumB.Gain.pval LumA.Loss.pval LumB.Loss.pval
LumA.Gain.pval <- c()
LumB.Gain.pval <- c()
LumA.Loss.pval <- c()
LumB.Loss.pval <- c()

pb = txtProgressBar(min = 0, max = nrow(cn.data.test), initial = 0, style = 3)
for(i in 1:length(row.names(cn.data.test))) {
  setTxtProgressBar(pb,i)
  # get data
  gene <- row.names(cn.data.test)[i]
  #print(i)
  #print(gene)
  her2e.dat <- as.numeric(her2e.data[gene,])
  luma.dat <- as.numeric(luma.data[gene,])
  lumb.dat <- as.numeric(lumb.data[gene,])
  
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


gene.test.df <- data.frame("gene"=row.names(cn.data.test), 
                           "LumA.Gain.pval"=LumA.Gain.pval,
                           "LumB.Gain.pval"=LumB.Gain.pval,
                           "LumA.Loss.pval"=LumA.Loss.pval,
                           "LumB.Loss.pval"=LumB.Loss.pval)

# adjust pval
gene.test.df$LumA.Gain.padj <- p.adjust(gene.test.df$LumA.Gain.pval, method = "fdr") # "bonferroni", "fdr"
gene.test.df$LumB.Gain.padj <- p.adjust(gene.test.df$LumB.Gain.pval, method = "fdr")
gene.test.df$LumA.Loss.padj <- p.adjust(gene.test.df$LumA.Loss.pval, method = "fdr")
gene.test.df$LumB.Loss.padj <- p.adjust(gene.test.df$LumB.Loss.pval, method = "fdr")
#View(gene.test.df)
# length(gene.test.df %>% filter(LumA.Gain.padj<=0.05) %>% pull(gene))
# length(gene.test.df %>% filter(LumB.Gain.padj<=0.05) %>% pull(gene))
# length(gene.test.df %>% filter(LumA.Loss.padj<=0.05) %>% pull(gene))
# length(gene.test.df %>% filter(LumB.Loss.padj<=0.05) %>% pull(gene))

save(gene.test.df, file= "data/METABRIC/4_CN/processed/Metabric_gene_cna_signif.RData")

# for the signfiicant genes, calc the difference in CNfreq and then visualize
genes <- loadRData("data/METABRIC/4_CN/processed/Metabric_gene_cna_signif.RData")
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

#cn.data.test[1:5,1:5]
gainloss.freqs <- cn.data.test %>% 
  rownames_to_column(var="gene") %>% 
  dplyr::select(gene)

# calc. loss/gain freqs per group
# her2e
gainloss.freqs$freq.loss.her2e <- apply(
  cn.data.test[,colnames(cn.data.test) %in% colnames(her2e.data)], 1, 
  function(x) {
    (length(which(x<=-1))/length(x))*-100
    }
  ) 

gainloss.freqs$freq.gain.her2e  <- apply(
  cn.data.test[,colnames(cn.data.test) %in% colnames(her2e.data)], 1, 
  function(x) {
    (length(which(x>=1))/length(x))*100
  }
)

# luma
gainloss.freqs$freq.loss.luma <- apply(
  cn.data.test[,colnames(cn.data.test) %in% colnames(luma.data)], 1, 
  function(x) {
    (length(which(x<=-1))/length(x))*-100
  }
) 

gainloss.freqs$freq.gain.luma <- apply(
  cn.data.test[,colnames(cn.data.test) %in% colnames(luma.data)], 1, 
  function(x) {
    (length(which(x>=1))/length(x))*100
  }
)

# lumb
gainloss.freqs$freq.loss.lumb <- apply(
  cn.data.test[,colnames(cn.data.test) %in% colnames(lumb.data)], 1, 
  function(x) {
    (length(which(x<=-1))/length(x))*-100
  }
) 

gainloss.freqs$freq.gain.lumb <- apply(
  cn.data.test[,colnames(cn.data.test) %in% colnames(lumb.data)], 1, 
  function(x) {
    (length(which(x>=1))/length(x))*100
  }
)

# save
save(gainloss.freqs,
     file = "data/METABRIC/4_CN/processed/CN_gainloss_frequencies.RData")

gene.freqs <- gainloss.freqs

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

print(p)
#boxplot(Freq.diff~PAM50,data=freq.diffs,col=c("#2176d5","#34c6eb"))
#p <- recordPlot()
#plot.new()
plot.list <- append(plot.list,list(p))



# save plots
pdf(file = plot.file, onefile = TRUE)#, height = 5, width = 5)

for(i in 1:length(plot.list)) { 
  print(i)
  print(plot.list[[i]])
}

dev.off()

# save text output
writeLines(txt.out, txt.file)
