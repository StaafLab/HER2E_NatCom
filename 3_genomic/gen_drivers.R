# Script: Driver mutations in the HER2E subtype (WGS based) - waterfall plot

#TODO: 

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "SCANB" # SCANB 

# set/create output directory for plots
output.path <- "output/plots/3_genomic/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/3_genomic/processed/",sep="")
dir.create(data.path)

# plot
plot.list <- list() # object to store plots; note: if the output is not in string format use capture.output()
plot.file <- paste(output.path,cohort,"_HER2n_driverWF.pdf",sep = "")

#packages
source("scripts/3_genomic/src/gen_functions.R")
source("scripts/4_CN/src/cn_functions.R")
library(ggplot2)
library(tidyverse)
library(readxl)
library(GenVisR)
library(reshape2)
#library(ggstatsplot)
#library(data.table)
library(grid)

#######################################################################
#######################################################################

# define list of driver genes/ cancer genes
#amp.drivers <- c("ERBB2", "CCND1", "ZNF703", "PAK1", "RPS6KB1", "MYC", "ZNF217", "MCL1", "MDM2", "PIK3CA", "CCNE1") # from TCGA

# use only basis paper - intersect with new data to see if all genes incl in this list
sub.indel.drivers <- as.data.frame(read_excel("data/BASIS/3_genomic/raw/Supplementary Table 14.Driver.Events.By.Mutation.Type.01052015.v2.xlsx", sheet = "Subs.indels")) %>% 
  pull(Gene) %>% 
  unique

amp.drivers <- as.data.frame(read_excel("data/BASIS/3_genomic/raw/Supplementary Table 14.Driver.Events.By.Mutation.Type.01052015.v2.xlsx", sheet = "CopyNumber")) %>% 
  mutate(Gene = if_else(Gene=="Chr8:(ZNF703/FGFR1)", "ZNF703", Gene)) %>% 
  pull(Gene) %>% 
  unique

driver.genes <- unique(c(sub.indel.drivers,amp.drivers))

#######################################################################

# ID key file
id.key <- as.data.frame(read_excel("./data/SCANB/3_genomic/raw/HER2_enriched_June23_ForJohan.xlsx", sheet = "Samples")) %>% 
  dplyr::select(c("Sample","Tumour")) %>% dplyr::rename(sample=Sample)

#######################################################################

# load driver data
indel.drivers <- as.data.frame(read_excel("./data/SCANB/3_genomic/raw/HER2_enriched_June23_ForJohan.xlsx", sheet = "AllCodingPindel")) %>% 
  mutate(VC = paste("indel_",VC, sep = "")) %>% 
  dplyr::rename(Tumour=Sample,gene=VD_Gene,variant_class=VC) %>% 
  left_join(id.key,by="Tumour") %>% relocate(sample,1) %>% 
  dplyr::select(c(sample,gene,variant_class)) %>% 
  filter(gene %in% driver.genes) %>% 
  distinct()

sub.drivers <- as.data.frame(read_excel("./data/SCANB/3_genomic/raw/HER2_enriched_June23_ForJohan.xlsx", sheet = "AllCodingASMD_CLPM")) %>% 
  mutate(VC = paste("sub_",VC, sep = "")) %>% 
  dplyr::rename(Tumour=Sample,gene=VD_Gene,variant_class=VC) %>% 
  left_join(id.key,by="Tumour") %>% relocate(sample,1) %>% 
  dplyr::select(c(sample,gene,variant_class)) %>% 
  filter(gene %in% driver.genes) %>% 
  distinct()

# # processed amplification data
# amp.dat <- loadRData("data/SCANB/4_CN/processed/CN_amp_genpos_genmap.RData")
# cn.drivers <- amp.dat[lengths(amp.dat$Gene_symbol) > 0,] %>% 
#   filter(!is.na(Gene_symbol)) %>% 
#   filter(as.character(Gene_symbol) %in% driver.genes) %>% 
#   filter(if_any(starts_with("S"), ~ . > 0))
# 
# # final matrix 
# cn.drivers.long <- data.frame()
# 
# # get it down to gene level
# for (sample in colnames(cn.drivers)[grepl("S",colnames(cn.drivers))]) {
#   # get sample data
#   sample.data <- cn.drivers[c(sample,"ProbeID","Gene_symbol")] #%>% filter(.[[sample]] > 0)
#   
#   # jump to next iteration if the sample has no driver amps
#   if (sort(unique(sample.data[[sample]]))[1]==0 & 
#       length(unique(sample.data[[sample]]))==1) { 
#     #print(paste("No amps in this sample: ",sample,sep=""))
#     next
#   }
#   #print(paste("AMP IN THIS SAMPLE: ",sample,sep=""))
#   for (gene in unique(sample.data$Gene_symbol)) {
#      # probes for a gene
#      gene.probes <- sample.data %>% 
#        filter(Gene_symbol == !!gene) %>% 
#        pull("ProbeID")
#      probe.count <- sample.data %>% 
#        filter(ProbeID %in% gene.probes) %>% 
#        filter(.[[sample]] > 0) %>% 
#        pull("ProbeID") %>% length(.)
#      if (probe.count >= length(gene.probes)/2) { # majority status called
#        cn.drivers.long <- rbind(cn.drivers.long,c(sample,gene,"amplified"))
#        
#      }
#    }
# }
# 
# # set col names
# names(cn.drivers.long) <- c("sample", "gene", "variant_class")
# 
# View(cn.drivers.long)
# HERE ADD NEW DATA replace cn.drivers.long
# amplification status
cna.genes <- loadRData("./data/SCANB/3_genomic/processed/ASCAT_genelevel.RData")
#View(cna.genes)
names(cna.genes) <- gsub("\\..*", "", names(cna.genes))
qc.scanb <- as.data.frame(read_excel("./data/SCANB/3_genomic/raw/HER2_enriched_June23_ForJohan.xlsx", sheet = "Samples"))
cna.genes <- cna.genes[names(cna.genes) %in% qc.scanb$Sample]
cna.driver.genes <- lapply(cna.genes, function(x) {
  #x <- cna.genes[[1]]
  filt.x <- x[x$gene %in% driver.genes & x$Amp==1,]
  filt.x <- filt.x[,c("sample","gene")]
  filt.x$sample <- gsub("\\..*", "", filt.x$sample)
  return(filt.x)
})
driv.amp <- do.call(rbind, cna.driver.genes)
rownames(driv.amp) <- NULL
colnames(driv.amp) <- c("sample","gene")
driv.amp$variant_class <- "amplified"

# make and save combined driver dataset for plotting
drivers.df <- do.call("rbind", list(driv.amp, sub.drivers, indel.drivers))
save(drivers.df, file = paste(data.path,"driver_mutations_all.RData",sep=""))
#load(paste(data.path,"driver_mutations_all.RData",sep=""))
#drivers.df[which(drivers.df$gene=="ERBB2"),]$sample # none mut with erbb2 high gex
################################################################################
# Waterfall plotting parameters
################################################################################

# plotting parameters
# 1. mainRecurCutoff accepts a numeric value between 0 and 1, and will only plot genes with mutations in x proportion of samples.
# 2. if there are specific genes of interest those can be specified directly via the plotGenes parameter. Input to plotGenes should be a character vector of a list of genes that are desireable to be shown and is case sensitive. 
# 3. plot only specific samples. This can be achieved via the parameter plotSamples
# 4. the maxGenes parameter will only plot the top x genes and takes an integer value. This is usefull for example if when using the mainRecurCutoff parameter a vector of genes have values at x cutoff and all of them are not desired. 
# 5. the rmvSilent parameter will remove all silent mutations from the data.

# combine
#mut.drivers <- rbind(sub.drivers,indel.drivers) #cn.drivers,
#mut.all <- rbind(sub.all,indel.all)

#colors <- c("#0e0421", "#d4136d", "#12e0dd", "#c70c0c", 
#                    "#2a18cc", "#0b9c32")
colors <- c("#8dd3c7",
            "#ffffb3",
            "#bebada",
            "#fb8072",
            "#80b1d3",
            "#fdb462",
            "#b3de69",
            "#fccde5")


################################################################################
# Waterfall plots
################################################################################

# datasets to be plotted
datasets <- list(drivers = drivers.df,
                 sub.drivers = sub.drivers,
                 indel.drivers = indel.drivers,
                 cn.drivers = driv.amp)

for (i in names(datasets)) { 
  
  # data
  data <- datasets[[i]] 
  mutation.priority <- as.character(unique(data$variant_class))
  custom.pallete <- colors[1:length(mutation.priority)]
  
  # title
  layer <- list(ggtitle(i))
  
  # plot # idea include all sample but only plot the 25 samples because toherwise the % mutatnt sidebar is not correct in relation to all 30 samples
  plot <- waterfall(data, 
                    fileType = "Custom", 
                    variant_class_order = mutation.priority,
                    mainGrid = TRUE,
                    mainPalette = custom.pallete,
                    main_geneLabSize = 15,
                    mainRecurCutoff = 0,
                    maxGenes = 10,
                    mainDropMut = TRUE, # drop unused mutation types from legend
                    #rmvSilent = TRUE,
                    out= "grob",
                    mutBurdenLayer = layer,
                    plotMutBurden = FALSE) #
  #plotSamples = c()
  
  #grid.draw(plot)
  
  # append to list
  plot.list <- append(plot.list,list(plot))
  
}

#######################################################################
#######################################################################

# save plots
pdf(file = plot.file, onefile = TRUE, height = 10, width = 15)

for (i in 1:length(plot.list)) {
  grid::grid.newpage()
  grid::grid.draw(plot.list[[i]])
  
  #print(plot.list[[i]])
}

dev.off()

