# Script: Checking PAM50 correlation

# TODO:

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "SCANB" # SCANB 

# set/create output directory for plots
output.path <- "output/plots/1_clinical/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/1_clinical/processed/",sep="")
dir.create(data.path)

#packages
#source("scripts/1_clinical/src/clin_functions.R")
library(ggplot2)
library(tidyverse)
library(reshape)
library(openxlsx)

# list to store plots
plot.list <- list()

#######################################################################
# Load data
#######################################################################

# load data
load("data/SCANB/1_clinical/raw/pam50ann_correlations_summarizedByMean_REL4.RData")
# load anno
clin.rel4 <- as.data.frame(read_excel("data/SCANB/1_clinical/raw/NPJ_release.xlsx"))

#######################################################################
# PAM50 centroid correlation in ERpHER2nHER2E samples
#######################################################################

# data processing already done: calc. mean (or median) for each sample for each pam50 group based on the 100 runs

# select subgroup data
anno <- clin.rel4 %>% 
    filter(Follow.up.cohort == TRUE) %>% 
    filter(ER=="Positive" & HER2=="Negative") %>% 
    filter(NCN.PAM50 == "Her2") 

data <- pam50ann_correlations %>% 
    filter(Assay %in% anno$GEX.assay) %>% 
    filter(majorityClass != "unclassified") %>% 
    filter(majorityNearest != majoritySecondBest) %>% 
    distinct() #S003571.l.r.m.c.lib.g.k2.a.t occurs 2 times for some reason, the data is identical (duplicated) so I exclude one

# creating modified dataframe
data_mod <- melt(data, measure.vars=colnames(data)[5:9])

# creating boxplot
plot.list <- append(plot.list,list(
    ggplot(data_mod) + 
    geom_boxplot(aes(x=variable, y=value,fill=as.factor(variable)),alpha=0.7, size=1.5, outlier.size = 5) +
    xlab("PAM50 Centroid") +
    ylab("Correlation") +
    ggtitle("PAM50 centroid correlation in ERpHER2nHER2E samples") +
    theme(plot.title = element_text(size = 25),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 25),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 25),
          legend.position = "none") +
    scale_fill_manual(values=c(meanLumA = "#2176d5", meanLumB = "#34c6eb", meanHer2 ="#d334eb", meanBasal ="#c41b0e", meanNormal="#64c25d"))))

#######################################################################
# Crosstable with 2nd best PAM50 match for all samples
#######################################################################

# select subgroup data
anno <- clin.rel4 %>% 
    filter(Follow.up.cohort == TRUE) %>% 
    filter(ER=="Positive" & HER2=="Negative") 

data <- pam50ann_correlations %>% 
    filter(Assay %in% anno$GEX.assay) %>% 
    filter(majorityClass != "unclassified") %>% 
    filter(majorityNearest != majoritySecondBest) %>% 
    distinct() #S003571.l.r.m.c.lib.g.k2.a.t occurs 2 times for some reason, the data is identical (duplicated) so I exclude one

confusion_matrix <- as.data.frame(table(data$majorityNearest,data$majoritySecondBest,dnn=c("majorityNearest","majoritySecondBest")))

plot.list <- append(plot.list,list(
    ggplot(data = confusion_matrix,
       mapping = aes(x = majorityNearest,
                     y = majoritySecondBest)) +
    scale_x_discrete(limits = rev(levels(confusion_matrix$majorityNearest))) +
    geom_tile(aes(fill = Freq)) +
    geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
    scale_fill_gradient(low = "#fee6ce",
                        high = "#e6550d",
                        trans = "log") +
    ggtitle("2nd best PAM50 match in ERpHER2n samples") +
    xlab("PAM50 Class") +
    ylab("2nd best PAM50 match") +
    theme(plot.title = element_text(size = 25),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 25),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 25),
          legend.position = "none")))

# small table for her2e
her2e_matrix <- confusion_matrix %>% 
    filter(majorityNearest == "Her2") %>% 
    droplevels() %>% 
    dplyr::select(-c(majorityNearest)) %>% 
    column_to_rownames(var="majoritySecondBest")
her2e_matrix <- as.data.frame(t(her2e_matrix)) %>% dplyr::select(-c(2))
write.xlsx(her2e_matrix, file = paste(data.path,"PAM50_HER2E_2ndbest.xlsx",sep=""),overwrite = TRUE)

#######################################################################
# Distinctiveness: Boxplot the difference between the best and second 
# best correlation matches
#######################################################################

# select subgroup data
anno <- clin.rel4 %>% 
    filter(Follow.up.cohort == TRUE) %>% 
    filter(ER=="Positive" & HER2=="Negative") 

data <- pam50ann_correlations %>% 
    filter(Assay %in% anno$GEX.assay) %>% 
    filter(majorityClass != "unclassified") %>% 
    filter(majorityNearest != majoritySecondBest) %>% 
    distinct() 

# dataframe best, diff to 2nd best
data <- data %>% rowwise() %>% mutate(Diff = abs(get(paste("mean",majorityNearest,sep=""))) - abs(get(paste("mean",majoritySecondBest,sep=""))))

# plot
plot.list <- append(plot.list,list(
    ggplot(data) + 
    geom_boxplot(aes(x=majorityNearest, y=Diff,fill=as.factor(majorityNearest)),alpha=0.7, size=1.5, outlier.size = 5) +
    xlab("PAM50 Class") +
    ylab("Correlation diff (1st-2nd)") +
    ggtitle("PAM50 class distinctiveness in ERpHER2n samples") +
    theme(plot.title = element_text(size = 25),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 25),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 25),
          legend.position = "none") +
    scale_fill_manual(values=c(LumA = "#2176d5", LumB = "#34c6eb", Her2 ="#d334eb", Basal ="#c41b0e", Normal="#64c25d"))))

#######################################################################
# HER2E specific distinctiveness: Boxplot for the her2e samples the difference between the best (her2e) and second best correlation matches
# so three groups in boxplot: basal, luma, lumb, normal
#######################################################################

# select subgroup data
anno <- clin.rel4 %>% 
    filter(Follow.up.cohort == TRUE) %>% 
    filter(ER=="Positive" & HER2=="Negative") %>% 
    filter(NCN.PAM50 == "Her2") 

data <- pam50ann_correlations %>% 
    filter(Assay %in% anno$GEX.assay) %>% 
    filter(majorityClass != "unclassified") %>% 
    filter(majorityNearest != majoritySecondBest) %>% 
    distinct() 

#table(data$majoritySecondBest)
# dataframe best, diff to 2nd best
data <- data %>% rowwise() %>% mutate(HER2Diff = abs(get(paste("mean",majorityNearest,sep=""))) - abs(get(paste("mean",majoritySecondBest,sep=""))))

# plot
plot.list <- append(plot.list,list(
    ggplot(data) + 
    geom_boxplot(aes(x=majoritySecondBest, y=HER2Diff,fill=as.factor(majoritySecondBest)),alpha=0.7, size=1.5, outlier.size = 5) +
    xlab("Second-best PAM50 Class") +
    ylab("Correlation diff (1st-2nd)") +
    ggtitle("PAM50-HER2E distinctiveness in ERpHER2nHER2E") +
    theme(plot.title = element_text(size = 25),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 25),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 25),
          legend.position = "none") +
    scale_fill_manual(values=c(LumA = "#2176d5", LumB = "#34c6eb", Her2 ="#d334eb", Basal ="#c41b0e", Normal="#64c25d"))))

# further investigate the luma,lumb,basal samples
# creating modified dataframe
data_mod <- data %>% dplyr::select(c(majoritySecondBest,HER2Diff,meanHer2))

# all samples
plot.list <- append(plot.list,list(
    ggplot(data_mod) + 
    geom_point(aes(x=meanHer2, y=HER2Diff,color=as.factor(majoritySecondBest)),alpha=0.7, size=5) +
    xlab("HER2E correlation (mean)") +
    ylab("Correlation diff (1st-2nd)") +
    ggtitle("HER2E distinctiveness: Subtype centroid correlations (ERpHER2nHER2E)") +
    theme(plot.title = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 25),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 25),
          legend.position = c(0.9,0.1),
          legend.text=element_text(size = 15),
          legend.title=element_text(size = 20)) +
    scale_color_manual("2nd Subtype",values=c(LumA = "#2176d5", LumB = "#34c6eb", Basal ="#c41b0e")) +
    scale_x_continuous(breaks=seq(0, 0.8, 0.1))))

# correlation between subtype and her2 centroid correlations
# luma
data_mod <- data %>% filter(majoritySecondBest == "LumA")

plot.list <- append(plot.list,list(
    ggplot(data_mod) + 
    geom_point(aes(x=meanHer2, y=meanLumA,color=as.factor(majoritySecondBest)),alpha=0.7, size=5) +
    xlab("HER2E correlation") +
    ylab("LumA correlation") +
    ggtitle("ERpHER2nHER2E: Luminal A and HER2E centroid correlations") +
    theme(plot.title = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 25),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 25),
          legend.position = "none") +
    scale_color_manual("Subtype",values=c(LumA = "#2176d5", LumB = "#34c6eb", Basal ="#c41b0e")) +
    scale_x_continuous(limits = c(0.2, 0.8),breaks=seq(0, 0.8, 0.1)) +
    scale_y_continuous(limits = c(0,0.6),breaks=seq(0, 0.6, 0.1))))

# lumb
data_mod <- data %>% filter(majoritySecondBest == "LumB")

plot.list <- append(plot.list,list(
    ggplot(data_mod) + 
    geom_point(aes(x=meanHer2, y=meanLumB,color=as.factor(majoritySecondBest)),alpha=0.7, size=5) +
    xlab("HER2E correlation") +
    ylab("LumB correlation") +
    ggtitle("ERpHER2nHER2E: Luminal B and HER2E centroid correlations") +
    theme(plot.title = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 25),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 25),
          legend.position = "none") +
    scale_color_manual("Subtype",values=c(LumA = "#2176d5", LumB = "#34c6eb", Basal ="#c41b0e")) +
    scale_x_continuous(limits = c(0.2, 0.8),breaks=seq(0, 0.8, 0.1)) +
    scale_y_continuous(limits = c(0,0.6),breaks=seq(0, 0.6, 0.1))))

# basal
data_mod <- data %>% filter(majoritySecondBest == "Basal")

plot.list <- append(plot.list,list(
    ggplot(data_mod) + 
    geom_point(aes(x=meanHer2, y=meanBasal,color=as.factor(majoritySecondBest)),alpha=0.7, size=5) +
    xlab("HER2E correlation") +
    ylab("Basal correlation") +
    ggtitle("ERpHER2nHER2E: Basal and HER2E centroid correlations") +
    theme(plot.title = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 25),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 25),
          legend.position = "none") +
    scale_color_manual("Subtype",values=c(LumA = "#2176d5", LumB = "#34c6eb", Basal ="#c41b0e")) +
    scale_x_continuous(limits = c(0.2, 0.8),breaks=seq(0, 0.8, 0.1)) +
    scale_y_continuous(limits = c(0,0.6),breaks=seq(0, 0.6, 0.1))))

# save plots
pdf(file = paste(output.path,cohort,"_PAM50correlation.pdf", sep=""), 
    onefile = TRUE, width = 21, height = 14.8) 

for (i in 1:length(plot.list)) {
    print(plot.list[[i]])
}

dev.off()
