# Script: Kataegis in SCANB and BASIS

#TODO: 

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "COMBINED"

# set/create output directory for plots
output.path <- "output/plots/4_CN/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/4_CN/processed/", sep="")
dir.create(data.path)
plot.file <- paste(output.path,cohort,"_HER2n_kat.pdf",sep = "")

txt.file <- paste(output.path,cohort,"_HER2n_kat_tests.txt", sep="")
txt.out <- c() # object to store text output

# store plots
plot.list <- list()

#packages
source("scripts/4_CN/src/cn_functions.R")
source("scripts/2_transcriptomic/src/tscr_functions.R")
library(ggplot2)
library(tidyverse)
library(readxl)

#######################################################################
#######################################################################

#kataegis data
kat.scanb <- loadRData("./data/SCANB/4_CN/raw/SCANB_ERp_collected_kataegis.RData")
kat.basis <- loadRData("./data/BASIS/4_CN/raw/BASIS_collected_kataegis.RData")

# get the annotation data for sample selection and comparison
scanb.anno <- loadRData("./data/SCANB/4_CN/processed/HER2Esample_ploidy.RData") %>% 
  mutate(PAM50="HER2E") %>% 
  dplyr::select(Sample,PAM50)
basis.anno <- loadRData("data/BASIS/1_clinical/raw/Summarized_Annotations_BASIS.RData") %>% 
  filter(ClinicalGroup == "ERposHER2neg" & PAM50_AIMS %in% c("LumA","LumB")) %>% 
  dplyr::select(c("sample_name","PAM50_AIMS")) %>% 
  dplyr::rename(Sample=sample_name,PAM50=PAM50_AIMS)


#######################################################################
#######################################################################

# get normal IDs in SCANB
id.key <- read_excel(
  "data/SCANB/3_genomic/raw/JVCimpression2022-11-09_ERpos_WGS_Batch1-3.xlsx") %>% 
  dplyr::select(c(SENT.TUMOR,SENT.TUMOR.aliquot,TUMOR.alias)) %>% 
  mutate(SENT.TUMOR.aliquot = gsub("\\.","_",SENT.TUMOR.aliquot))
id.df <- sapply(list(names(kat.scanb)), function(x) {sub("_vs_.*", "", x)}) %>% 
  as.data.frame() %>% 
  dplyr::rename(Old.id=1)
# 1. convert epb IDs to normal sample IDs
id.df$New.id <- id.key$SENT.TUMOR[match(
  id.df$Old.id, id.key$SENT.TUMOR.aliquot)] 
# 2. convert the other IDs to normal sample IDs
id.df$New.id.2 <- id.key$SENT.TUMOR[match(
  id.df$Old.id, id.key$TUMOR.alias)] # 2. convert the other IDs to normal sample IDs
id.df$New.id <- ifelse(is.na(id.df$New.id), id.df$New.id.2, id.df$New.id)
id.df$New.id.2 <- NULL
# name
names(kat.scanb) <- id.df$New.id

# get normal IDs in BASIS
names(kat.basis) <- sapply(list(names(kat.basis)), 
                           function(x) {sub("PD(.*?)[A-Za-z].*", "PD\\1", x)})

#######################################################################

# get my sample data
kat.scanb <- kat.scanb[scanb.anno$Sample]
kat.basis <- kat.basis[basis.anno$Sample]

# get kataegis counts per tumor 
kat.scanb.counts <- lapply(kat.scanb, function(x) {
  x <- x[which(x$confidence>=1),] # filter by confidence
  return(nrow(x)) # number of events
})
# Convert the list to a data frame
kat.scanb.counts <- data.frame(
  SampleID = names(kat.scanb.counts),
  Kat_events = unlist(kat.scanb.counts)
) %>% 
  mutate(PAM50 = "HER2E")

# get kataegis counts per tumor 
kat.basis.counts <- lapply(kat.basis, function(x) {
  x <- x[which(x$confidence>=1),] # filter by confidence
  return(nrow(x)) # number of events
})
# Convert the list to a data frame
kat.basis.counts <- data.frame(
  SampleID = names(kat.basis.counts),
  Kat_events = unlist(kat.basis.counts)
)
kat.basis.counts$PAM50 <- basis.anno$PAM50[match(kat.basis.counts$SampleID,basis.anno$Sample)]

kat.counts <- rbind(kat.basis.counts,kat.scanb.counts)

# binary column
kat.counts <- kat.counts %>% mutate(Kat_binary = ifelse(Kat_events == 0, 0, 1))
save(kat.counts,file = "./data/SCANB/4_CN/processed/CN_kataegis.RData")
#######################################################################
#######################################################################

p <- ggplot(kat.counts, aes(x=PAM50,y=as.numeric(Kat_events),fill=PAM50)) +
  geom_boxplot(size=2.5, outlier.size = 7) +
  ylab("Kataegis events") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",linewidth=2),
        axis.ticks = element_line(colour = "black", linewidth = 2)) +
  scale_fill_manual(values=setNames(c("#d334eb","#2176d5","#34c6eb"),
                                    c("HER2E","LumA","LumB"))) +  
  scale_y_continuous(breaks = scales::breaks_pretty(10))


#print(p)
plot.list <- append(plot.list,list(p))#

p <- ggplot(kat.counts, aes(x=PAM50,y=as.numeric(Kat_events),fill=PAM50)) +
  geom_boxplot(size=2.5, outlier.size = 7) +
  ylab("Kataegis events") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",linewidth=2),
        axis.ticks = element_line(colour = "black", linewidth = 2)) +
  scale_fill_manual(values=setNames(c("#d334eb","#2176d5","#34c6eb"),
                                    c("HER2E","LumA","LumB"))) +  
  scale_y_continuous(limits = c(0,25),breaks = scales::breaks_pretty(10))


#print(p)
plot.list <- append(plot.list,list(p))#
#######################################################################
# stats

# wilcox
txt.out <- append(txt.out,c("Testing Kataegis event counts"))
res <- mwu_test(kat.counts[which(kat.counts$PAM50 == "HER2E"),]$Kat_events,
         kat.counts[which(kat.counts$PAM50 == "LumA"),]$Kat_events)
txt.out <- append(txt.out,c(capture.output(res)))
res <- mwu_test(kat.counts[which(kat.counts$PAM50 == "HER2E"),]$Kat_events,
         kat.counts[which(kat.counts$PAM50 == "LumB"),]$Kat_events)
txt.out <- append(txt.out,c(capture.output(res)))

# freqs
her2e.dat <- kat.counts[kat.counts$PAM50 == "HER2E",]$Kat_binary
her2e.freq <- (length(her2e.dat[her2e.dat != 0])/length(her2e.dat))*100
luma.dat <- kat.counts[kat.counts$PAM50 == "LumA",]$Kat_binary
luma.freq <- (length(luma.dat[luma.dat != 0])/length(luma.dat))*100
lumb.dat <- kat.counts[kat.counts$PAM50 == "LumB",]$Kat_binary
lumb.freq <- (length(lumb.dat[lumb.dat != 0])/length(lumb.dat))*100

# chi2
txt.out <- append(txt.out,c("Testing Kataegis binary"))
txt.out <- append(txt.out,c("Kataegis frequencies:",
                            "  HER2E =",her2e.freq,
                            "  LumA =",luma.freq,
                            "  LumB =",lumb.freq))

res <- fisher.test(table(
  kat.counts[which(
    kat.counts$PAM50 %in% c("HER2E","LumA")),c("PAM50", "Kat_binary")]))
txt.out <- append(txt.out,c(capture.output(res)))
res <- fisher.test(table(
  kat.counts[which(
    kat.counts$PAM50 %in% c("HER2E","LumB")),c("PAM50", "Kat_binary")]))
txt.out <- append(txt.out,c(capture.output(res)))

#######################################################################
#######################################################################

# save plots
pdf(file = plot.file, onefile = TRUE)#, height = 5, width = 5)

for(i in 1:length(plot.list)) { 
  print(i)
  print(plot.list[[i]])
}

dev.off()

# save text output
writeLines(txt.out, txt.file)
