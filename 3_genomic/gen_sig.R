# Script: Mutational and rearrangement signatures in the HER2E subtype (WGS based)

# TODO: still need to fix the statsitcs part
# - 

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run (for plot naming)
cohort <- "COMBINED" 

# set/create output directory for plots
output.path <- "output/plots/3_genomic/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/3_genomic/processed/",sep="")
dir.create(data.path)

# plot
plot.list <- list() # object to store plots; note: if the output is not in string format use capture.output()
plot.file <- paste(output.path,cohort,"_WGS_signatures.pdf",sep = "")

txt.out <- c() # object to store text output
txt.file <- paste(output.path,cohort,"_WGS_signatures.txt", sep="")

#packages
#source("scripts/3_genomic/src/mut_functions.R")
source("scripts/4_CN/src/cn_functions.R")
library(ggplot2)
library(tidyverse)
library(readxl)

pdf(file = plot.file, onefile = TRUE, height = 5, width = 5)

################################################################################
# Mutational Signatures
################################################################################

# load the raw signature counts
mut.scanb <- as.data.frame(read_excel("./data/SCANB/3_genomic/raw/HER2_enriched_June23_ForJohan.xlsx", sheet = "SBSsig")) %>% dplyr::select(-c("Sample_pair","...3","...15"))
names(mut.scanb) <- gsub('SBS', 'S', colnames(mut.scanb))
mut.scanb <- column_to_rownames(mut.scanb, var = "Sample")

# convert SCANB to proportion per sample
mut.scanb <- as.data.frame(t(apply(t(mut.scanb),2,function(x) (x/sum(x, na.rm=TRUE))))) 

# load the raw signature counts
mut.basis <- loadRData("data/BASIS/1_clinical/raw/Summarized_Annotations_BASIS.RData")
# I have to define 2 groups: ERpHER2n LumA; ERpHER2n LumA
mut.basis <- as.data.frame(mut.basis) %>% 
    dplyr::select(c(sample_name,final.ER,final.HER2,
             ClinicalGroup,PAM50_AIMS,MutSigProportions)) %>%
    mutate(Subtype = 
               case_when(final.ER=="positive" & final.HER2=="negative" & PAM50_AIMS=="LumA" ~ "LumA",
                         final.ER=="positive" & final.HER2=="negative" & PAM50_AIMS=="LumB" ~ "LumB")) %>% 
    dplyr::filter(Subtype %in% c("LumA","LumB")) %>% 
    dplyr::select(Subtype,MutSigProportions)

mut.basis$MutSigProportions <- as.data.frame(mut.basis$MutSigProportions)
#table(mut.basis$Subtype)

# check which signatures are present in both cohorts
common.mut.sigs <- intersect(colnames(mut.basis$MutSigProportions),colnames(mut.scanb))

# exclude non-shared ones
mut.scanb <- mut.scanb[,which(names(mut.scanb) %in% common.mut.sigs)]
mut.basis$MutSigProportions <- mut.basis$MutSigProportions[,which(
  names(mut.basis$MutSigProportions) %in% common.mut.sigs)]

# add subtype info
mut.scanb <- mut.scanb %>% mutate(Subtype = "HER2E")

# get basis data in right format
mut.basis <- as.data.frame(do.call(cbind, mut.basis))
names(mut.basis) <- gsub("MutSigProportions.","", names(mut.basis))

# put together in an object ready for plotting 
mut.data <- rbind(mut.basis,mut.scanb)

################################################################################
# plot 

# luminal subtypes
group.colors <- c(LumA = "#0f1bbd", LumB = "#09d3e6", HER2E ="#b609e6")
plot.data <- mut.data  #%>% filter(Subtype != "ERpHER2p")

for(i in c(1:7)) {
  
    plot <- ggplot(plot.data, aes(
      x=as.factor(Subtype),y=plot.data[,i+1],fill=as.factor(Subtype))) +
                geom_boxplot(alpha=0.7, size=1.5, outlier.size = 5) +
                xlab("Subtype") +
                ylab("Proportion") +
                ylim(c(0,1)) +
                ggtitle(paste("Mutational signature ",colnames(plot.data)[i+1],sep="")) +
                theme_bw() +
                theme(aspect.ratio=1/1,
                      legend.position = "none",
                      panel.border = element_blank(), 
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black",linewidth = 2),
                      axis.ticks = element_line(colour = "black", linewidth = 1),
                      axis.ticks.length=unit(0.2, "cm")) +
                scale_fill_manual(values=group.colors)
    
    plot.list <- append(plot.list,list(plot))
    print(plot)
}

################################################################################
# Statistics
# vars
mut.data$Subtype <- as.factor(mut.data$Subtype)
for (signature in colnames(mut.data)[-1]) {
  her2e.data <- subset(mut.data, Subtype == "HER2E", select=c(signature))
  for (subtype in c("LumA","LumB")) {
    # signature data for subtype
    comp.data <- subset(mut.data, Subtype == subtype, select=c(signature))
    
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
    # save result
    txt.out <- append(txt.out,c(signature,subtype,capture.output(res)))
  }
}

  
################################################################################
# Rearrangement Signatures
################################################################################

# load the raw signature counts
rearr.scanb <- as.data.frame(read_excel("./data/SCANB/3_genomic/raw/HER2_enriched_June23_ForJohan.xlsx", sheet = "RearrangementSig")) %>% dplyr::select(-c("Sample_pair","Normal","HER2enriched"))
names(rearr.scanb) <- gsub('RefSigR', 'RS', colnames(rearr.scanb))
rearr.scanb <- column_to_rownames(rearr.scanb, var = "Sample")

# convert SCANB to proportion per sample
rearr.scanb <- as.data.frame(t(apply(t(rearr.scanb),2,function(x) (x/sum(x, na.rm=TRUE))))) 

# load the raw signature counts
rearr.basis <- loadRData("data/BASIS/1_clinical/raw/Summarized_Annotations_BASIS.RData")
rearr.basis <- rearr.basis %>% 
    dplyr::select(c(sample_name,final.ER,final.HER2,
             ClinicalGroup,PAM50_AIMS,RSproportions)) %>%
    mutate(Subtype = 
               case_when(final.ER=="positive" & final.HER2=="negative" & PAM50_AIMS=="LumA" ~ "LumA",
                         final.ER=="positive" & final.HER2=="negative" & PAM50_AIMS=="LumB" ~ "LumB")) %>% 
    filter(Subtype %in% c("LumA","LumB")) %>% 
    dplyr::select(Subtype,RSproportions)
rearr.basis$RSproportions <- as.data.frame(rearr.basis$RSproportions)

# pool= RS6 = RS6a + RS6b
rearr.scanb <- rearr.scanb %>% mutate(RS6 = RS6a+RS6b) %>% dplyr::select(-c(RS6a,RS6b))
#str(rearr.scanb)

# check which signatures are present in both cohorts
common.rearr.sigs <- intersect(colnames(rearr.basis$RSproportions),colnames(rearr.scanb))

# exclude non-shared ones
rearr.scanb <- rearr.scanb[,which(names(rearr.scanb) %in% common.rearr.sigs)]
rearr.basis$RSproportions <- rearr.basis$RSproportions[,which(
  names(rearr.basis$RSproportions) %in% common.rearr.sigs)]

# add subtype
rearr.scanb <- rearr.scanb %>% 
    mutate(Subtype = "HER2E")

# get basis data in right format
rearr.basis <- as.data.frame(do.call(cbind, rearr.basis))
names(rearr.basis) <- gsub("RSproportions.","", names(rearr.basis))

# put together in an object ready for plotting 
rearr.data <- rbind(rearr.basis,rearr.scanb)

################################################################################
# plot
group.colors <- c(LumA = "#0f1bbd", LumB = "#09d3e6", HER2E ="#b609e6")

for(i in c(1:6)) {
    plot <- ggplot(rearr.data, aes(x=as.factor(Subtype),y=rearr.data[,i+1],fill=as.factor(Subtype))) +
        geom_boxplot(alpha=0.7, size=1.5, outlier.size = 5) +
        xlab("Subtype") +
        ylab("Proportion") +
        ylim(c(0,1)) +
        ggtitle(paste("Rearrangement signature ",colnames(rearr.data)[i+1],sep="")) +
      theme_bw() +
      theme(aspect.ratio=1/1,
            legend.position = "none",
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black",linewidth = 2),
            axis.ticks = element_line(colour = "black", linewidth = 1),
            axis.ticks.length=unit(0.2, "cm")) +
        scale_fill_manual(values=group.colors)
    
    plot.list <- append(plot.list,list(plot))
    print(plot)
}

################################################################################

# Statistics
# vars
rearr.data$Subtype <- as.factor(rearr.data$Subtype)
for (signature in colnames(rearr.data)[-1]) {
  her2e.data <- subset(rearr.data, Subtype == "HER2E", select=c(signature))
  for (subtype in c("LumA","LumB")) {
    # signature data for subtype
    comp.data <- subset(rearr.data, Subtype == subtype, select=c(signature))
    
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
    # save result
    txt.out <- append(txt.out,c(signature,subtype,capture.output(res)))
  }
}

#######################################################################
#######################################################################

# save plots
# pdf(file = plot.file, onefile = TRUE, height = 5, width = 5)

# for(i in 1:length(plot.list)) { 
#   print(i)
#   print(plot.list[[i]])
# }

dev.off()

# save text output
writeLines(txt.out, txt.file)
