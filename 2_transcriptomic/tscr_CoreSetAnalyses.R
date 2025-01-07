# Script: 
# Core gene set analyses

#TODO: 

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

#-------------------------------------------------------------------------------
# packages
source("scripts/4_CN/src/cn_functions.R")
source("./scripts/2_transcriptomic/src/tscr_functions.R")
library(ggplot2)
library(tidyverse)
library(reshape2)
#library(data.table)
library(purrr)
library(readxl)
library(IRanges)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(Repitools)
library(openxlsx)
library(enrichR)
library(janitor)
library(ggvenn)
library(GenVisR)
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
DE.res.scanb <- "data/SCANB/2_transcriptomic/processed/DE_results.RData"
DE.res.metabric <- "data/METABRIC/2_transcriptomic/processed/DE_results.RData"
CS.cluster1 <- "data/SCANB/2_transcriptomic/processed/DEGs_HER2E_Sub1_network_nCon4_corCut0.5_GENES.txt"
CS.cluster2 <-"data/SCANB/2_transcriptomic/processed/DEGs_HER2E_Sub2_network_nCon4_corCut0.5_GENES.txt"
TFBS.transfak.cluster1 <- "data/SCANB/2_transcriptomic/processed/TF_Sub1_Transfac_tfsearch_TFBS_significant_analysis_file_100000_p(1.0).txt"
TFBS.transfak.cluster2 <- "data/SCANB/2_transcriptomic/processed/TF_Sub2_Transfac_tfsearch_TFBS_significant_analysis_file_100000_p(1.0).txt"
TFBS.transfak.all <- "data/SCANB/2_transcriptomic/processed/TF_Lenart_Transfac_tfsearch_TFBS_significant_analysis_file_100000_p(1.0).txt"
TFBS.jaspar.all <- "data/SCANB/2_transcriptomic/processed/TF_Lenart_JASPAR_new_tfsearch_TFBS_significant_analysis_file_100000_p(1.0).txt"
scanb.anno <- "./data/SCANB/1_clinical/processed/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data_ERpHER2n.RData"
scanb.gex <- "data/SCANB/2_transcriptomic/raw/genematrix_noNeg.Rdata"
gene.anno.path <- "data/SCANB/1_clinical/raw/Gene.ID.ann.Rdata"


# output
plot.file <- "output/plots/4_CN/COMBINED_coreDEGs.pdf"
txt.file <- "output/plots/4_CN/COMBINED_coreDEGs.txt"

# Open pdf file
pdf(file= "output/plots/4_CN/COMBINED_coreDEGs.pdf",width = 8, height = 8)

# create a 2X2 grid
par(mfrow= c(2,2))

#-------------------------------------------------------------------------------
# storing objects 
plot.list <- list() # object to store plots; note: if the output is not in string format use capture.output()
txt.out <- c() # object to store text output
# snippet to save plots (paste to end of script)
# pdf(file = plot.file, onefile = TRUE, height = 5, width = 5)
# for(i in 1:length(plot.list)) { 
#   print(i)
#   print(plot.list[[i]])
# }
# dev.off()
# snippet to save text output (paste to end of script)
#writeLines(txt.out, txt.file)
#-------------------------------------------------------------------------------
# snipped to take exection time
#start.time <- Sys.time()
# CODE
#end.time <- Sys.time()
#time.taken <- end.time - start.time
#time.taken

################################################################################
# required data
################################################################################

# load data
DE.res.scanb <- loadRData(DE.res.scanb)
DE.res.metabric <- loadRData(DE.res.metabric)

################################################################################
# Core DEGs definition and visualization
################################################################################

pam50.genes <- c("ACTR3B", "ANLN", "BAG1", "BCL2", "BIRC5", "BLVRA", "CCNB1", "CCNE1", "CDC20", "CDC6", "CDH3", "CENPF", "CEP55", "CXXC5", "EGFR", "ERBB2", "ESR1", "EXO1", "FGFR4", "FOXA1", "FOXC1", "GPR160", "GRB7", "KIF2C", "KRT14", "KRT17", "KRT5", "MAPT", "MDM2", "MELK", "MIA", "MKI67", "MLPH", "MMP11", "MYBL2", "MYC", "NAT1", "NDC80", "NUF2", "ORC6L", "PGR", "PHGDH", "PTTG1", "RRM2", "SFRP1", "SLC39A6", "TMEM45B", "TYMS", "UBE2C", "UBE2T")

# define DEG set
DEGs.scanb <- DE.res.scanb %>% 
  filter(Her2.LumA.padj <= 0.05) %>% 
  filter(Her2.LumB.padj <= 0.05) %>% 
  rownames_to_column("Gene") %>% pull(Gene)

DEGs.metabric <- DE.res.metabric %>% 
  filter(Her2.LumA.padj <= 0.05) %>% 
  filter(Her2.LumB.padj <= 0.05) %>% 
  rownames_to_column("Gene") %>% pull(Gene) 

# core gex
DEGs.core <- intersect(DEGs.scanb,DEGs.metabric)
pam50.DEGs.core <- intersect(pam50.genes,DEGs.core)

dat <- list("SCAN-B DEGs"=DEGs.scanb,"Metabric DEGs"=DEGs.metabric)
plot <- ggvenn(dat, 
  fill_color = c("#EFC000FF", "#CD534CFF"),
  stroke_size = 1, set_name_size = 5,
  show_percentage = FALSE
)
#grid::grid.newpage()
grid::grid.draw(plot)

# define DEG set
DEGs.scanb.la <- DE.res.scanb %>% 
  filter(Her2.LumA.padj <= 0.05) %>% 
  rownames_to_column("Gene") %>% pull(Gene)
DEGs.scanb.lb <- DE.res.scanb %>% 
  filter(Her2.LumB.padj <= 0.05) %>% 
  rownames_to_column("Gene") %>% pull(Gene)

DEGs.metabric.la <- DE.res.metabric %>% 
  filter(Her2.LumA.padj <= 0.05) %>% 
  rownames_to_column("Gene") %>% pull(Gene) 

DEGs.metabric.lb <- DE.res.metabric %>% 
  filter(Her2.LumB.padj <= 0.05) %>% 
  rownames_to_column("Gene") %>% pull(Gene) 

dat <- list("SCAN-B LumA-DEGs"=DEGs.scanb.la,"SCAN-B Lumb-DEGs"=DEGs.scanb.lb, "Metabric LumA-DEGs"=DEGs.metabric.la,"Metabric LumB-DEGs"=DEGs.metabric.lb)

plot <- ggvenn(dat, 
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
       stroke_size = 1, set_name_size = 5,
       show_percentage = FALSE
)
#grid::grid.newpage()
grid::grid.draw(plot)
################################################################################
# Gene set enrichment of core set
################################################################################
# set database parameters
setEnrichrSite("Enrichr") 
dbs <- listEnrichrDbs() # see avaiable dbs and select 

dbs <- c("WikiPathway_2023_Human")#"GO_Biological_Process_2023")#, "WikiPathway_2021_Human", "KEGG_2021_Human", "GO_Molecular_Function_2023")
res.core <- as.data.frame(enrichr(DEGs.core, dbs)[[1]])  %>% 
  mutate(Comp="All") %>% 
  mutate(Gene_count = as.numeric(sapply(strsplit(Overlap, "/"), "[[", 1))) %>% 
  arrange(dplyr::desc(Gene_count),Adjusted.P.value) %>% 
  mutate(Adjusted.P.value= round(Adjusted.P.value,3)) %>% 
  dplyr::slice(1:3)

################################################################################
# Cytoscape clusters - GSEA and HER2-low check
################################################################################

CS.clust.1 <- read.table(CS.cluster1) %>% pull(V1)
CS.clust.2 <- read.table(CS.cluster2) %>% pull(V1)
#length(CS.clust.1)+length(CS.clust.2)
# 1. GSEA
res.CS.clust.1 <- as.data.frame(enrichr(CS.clust.1, dbs)[[1]])
save(res.CS.clust.1, file = "./output/supplementary_data/HER2E_cluster1_PWE_WP2023.RData")
res.CS.clust.1 <- res.CS.clust.1  %>% 
  mutate(Comp="Sub1") %>% 
  mutate(Gene_count = as.numeric(sapply(strsplit(Overlap, "/"), "[[", 1))) %>% 
  arrange(dplyr::desc(Gene_count),Adjusted.P.value) %>% 
  mutate(Adjusted.P.value= round(Adjusted.P.value,3)) %>% 
  dplyr::slice(1:3)

res.CS.clust.2 <- as.data.frame(enrichr(CS.clust.2, dbs)[[1]])
save(res.CS.clust.2, file = "./output/supplementary_data/HER2E_cluster2_PWE_WP2023.RData")

res.CS.clust.2 <- res.CS.clust.2  %>% 
  mutate(Comp="Sub2") %>% 
  mutate(Gene_count = as.numeric(sapply(strsplit(Overlap, "/"), "[[", 1))) %>% 
  arrange(dplyr::desc(Gene_count),Adjusted.P.value) %>% 
  mutate(Adjusted.P.value= round(Adjusted.P.value,3)) %>% 
  dplyr::slice(1:3)

res.all <- rbind(res.core,res.CS.clust.1,res.CS.clust.2) %>% 
  mutate(Term = gsub("\\s*\\([^\\)]+\\)","",as.character(Term)))
  
res.all$Comp <- as.factor(res.all$Comp)
res.all$Term <- strtrim(res.all$Term, 52)

# labeller function for plot
deg.sets <- list("Core" = paste("Core DEGs (n=",length(DEGs.core),")",sep=""),
                 "Sub1" = paste("Sub1 (n=",length(CS.clust.1),")",sep=""),
                 "Sub2" = paste("Sub2 (n=",length(CS.clust.2),")",sep=""))
set_labeller <- function(variable,value){
  return(deg.sets[value])
}
# plot enrichment analysis
pwplot <- ggplot(res.all, aes(y=Gene_count,
                  x=Term,
                  fill=Adjusted.P.value)) +
  geom_bar(stat="identity") +
  coord_flip() +
  theme_bw() +
  labs(fill = "Adj. p-value") +
  scale_y_continuous(limits = c(0,65),breaks = seq(0,65,5)) +
  scale_fill_gradientn(limits=c(0,0.5), colors=c("red","blue")) +
  facet_grid(rows = vars(Comp),drop=TRUE,scales="free",space="free",labeller=set_labeller)

#grid::grid.newpage()
grid::grid.draw(pwplot)



# 2. HER2-low check in ERpHER2n stratified by PAM50
#  to understand if these genes "target" HER2-low cancers as well. 

# metagene-like expression (mean FPKM) for each pam50 subtype -> h2low yes no
# load annotation data and select subgroup data
anno <- loadRData(file=scanb.anno) %>%
  filter(Follow.up.cohort==TRUE) %>% 
  dplyr::rename(sampleID = GEX.assay, PAM50 = NCN.PAM50) %>% 
  dplyr::select(sampleID, PAM50, HER2_Low)

# load gene anno data to convert IDs
gene.anno <- loadRData(gene.anno.path) %>% as.data.frame() 

# load gex data
gex.data.core <- loadRData(scanb.gex) %>% as.data.frame() %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  mutate(geneID = gene.anno[which(gene.anno$Gene.ID==ensembl_gene_id),][["HGNC"]]) %>% 
  dplyr::select(-c(ensembl_gene_id)) %>% 
  drop_na(geneID) %>% 
  distinct(geneID,.keep_all = TRUE) %>% 
  column_to_rownames("geneID") %>%
  dplyr::select(any_of(anno$sampleID)) %>% 
  filter(row.names(.) %in% DEGs.core) %>% 
  select_if(~ !any(is.na(.))) 

# log transform FPKM data
gex.data.core <- as.data.frame(log2(gex.data.core + 1))
# z-transform
gex.data.core <- as.data.frame(t(apply(gex.data.core, 1, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y))))) 
gex.data.clust1 <- gex.data.core %>% 
  filter(row.names(.) %in% CS.clust.1)

gex.data.clust2 <- gex.data.core %>% 
  filter(row.names(.) %in% CS.clust.2)

gexdat.list <- list("Core_set" = gex.data.core, "Cluster_1" = gex.data.clust1,
                    "Cluster_2" = gex.data.clust2)
res.list <- list()
for(i in 1:length(gexdat.list)) {
  # gene set
  set <- names(gexdat.list)[i]
  set.data <- gexdat.list[[i]]
  # calc. the score for each sample
  result <- as.data.frame(apply(set.data, 2, mean)) %>% 
    dplyr::rename(!!set := "apply(set.data, 2, mean)") %>% 
    rownames_to_column(var="sampleID")
  res.list <- append(res.list,result)
}
res.df <- do.call(cbind,res.list)
res.df <- res.df[,c(1,2,4,6)] %>% as.data.frame %>% 
  dplyr::mutate_at(c("Core_set","Cluster_1","Cluster_2"), as.numeric)
# add pam50 anno
res.df$PAM50 <- anno$PAM50[match(res.df$sampleID,anno$sampleID)]
res.df$HER2_Low <- anno$HER2_Low[match(res.df$sampleID,anno$sampleID)]
res.df$HER2_Low <- ifelse(res.df$HER2_Low==1, "HER2Low_yes","HER2Low_no" )

# plot boxplots
# general 

boxplot(res.df$Core_set ~ res.df$PAM50, ylab="mean FPKM", main="Core set expression")
        #col = c("#d334eb","#2176d5","#34c6eb"))
boxplot(res.df$Cluster_1 ~ res.df$PAM50, ylab="mean FPKM", main="Sub1 cluster expression")
boxplot(res.df$Cluster_2 ~ res.df$PAM50, ylab="mean FPKM", main="Sub2 cluster expression")
res.df$HER2_Low[res.df$HER2_Low=="HER2Low_yes"] = "H2Low.1"
res.df$HER2_Low[res.df$HER2_Low=="HER2Low_no"] = "H2Low.0"

boxplot(res.df$Core_set ~ interaction(res.df$HER2_Low,res.df$PAM50), 
        ylab="mean FPKM", main="Core set expression vs. HER2-low",
        cex.main=1,las=2,xlab=NULL)
boxplot(res.df$Cluster_1 ~ interaction(res.df$HER2_Low,res.df$PAM50), 
        ylab="mean FPKM", main="Sub1 cluster expression vs. HER2-low",
        cex.main=1,las=2,xlab=NULL)
boxplot(res.df$Cluster_2 ~ interaction(res.df$HER2_Low,res.df$PAM50), 
        ylab="mean FPKM", main="Sub2 cluster expression vs. HER2-low",
        cex.main=1,las=2,xlab=NULL)

################################################################################
# TFBS analysis results - by Sunny
################################################################################

TFBS.transfak.cluster1 <- read.table(TFBS.transfak.cluster1,header=TRUE,sep = "\t")
head(TFBS.transfak.cluster1)
TFBS.transfak.cluster2 <- read.table(TFBS.transfak.cluster2,header=TRUE,sep = "\t")
head(TFBS.transfak.cluster2)
TFBS.transfak.all <- read.table(TFBS.transfak.all,header=TRUE,sep = "\t")
head(TFBS.transfak.all)
TFBS.jaspar.all <- read.table(TFBS.jaspar.all,header=TRUE,sep = "\t")
head(TFBS.jaspar.all)
#View(TFBS.transfak.all)
################################################################################
# Expression of top TFs
################################################################################

top.TFs.list <- list("TFBS.transfak.all"=TFBS.transfak.all[["TF"]][1:10],
                "TFBS.jaspar.all"=TFBS.jaspar.all[["TF"]][1:10],
                "TFBS.transfak.cluster1"=TFBS.transfak.cluster1[["TF"]][1:10],
                "TFBS.transfak.cluster2"=TFBS.transfak.cluster2[["TF"]][1:10])

# plot expr of top 10 TFs (her2e vs. luma vs. lumb)
top.TFs <- unique(unname(unlist(top.TFs.list)))

# get gex data
tf.gex.dat <- loadRData(scanb.gex) %>% as.data.frame() %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  mutate(geneID = gene.anno[which(gene.anno$Gene.ID==ensembl_gene_id),][["HGNC"]]) %>% 
  dplyr::select(-c(ensembl_gene_id)) %>% 
  drop_na(geneID) %>% 
  distinct(geneID,.keep_all = TRUE) %>% 
  column_to_rownames("geneID") %>%
  dplyr::select(any_of(anno$sampleID)) %>% 
  filter(row.names(.) %in% top.TFs) %>% 
  select_if(~ !any(is.na(.))) 
#View(tf.gex.dat)
# log transform FPKM data
tf.gex.dat <- as.data.frame(log2(tf.gex.dat + 1))
# z-transform
tf.gex.dat <- as.data.frame(t(apply(tf.gex.dat, 1, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y))))) 

for (i in 1:nrow(tf.gex.dat)) {
  tf <- row.names(tf.gex.dat)[i]
  dat <- t(tf.gex.dat[tf,]) %>% as.data.frame() %>% rownames_to_column(var="sampleID")
  dat$PAM50 <- anno$PAM50[match(dat$sampleID, anno$sampleID)]

  boxplot(dat[[tf]] ~ dat$PAM50,ylab=tf,
          # returns true if string is in list
          main=paste(names(top.TFs.list)[sapply(top.TFs.list, function(x) { 
    return(is.element(tf, x))})], collapse = ' & '))
}

################################################################################
# Waterfall plot of top TFs (mut status and amp)
################################################################################

# get mut and amp status for all available tf genes
top.TFs <- unique(unname(unlist(top.TFs.list)))

id.key <- as.data.frame(read_excel("./data/SCANB/3_genomic/raw/HER2_enriched_June23_ForJohan.xlsx", sheet = "Samples")) %>% 
  dplyr::select(c("Sample","Tumour")) %>% dplyr::rename(sample=Sample)

# sample; gene; mutation
indel.drivers <- as.data.frame(read_excel("./data/SCANB/3_genomic/raw/HER2_enriched_June23_ForJohan.xlsx", sheet = "AllCodingPindel")) %>% 
  mutate(VC = paste("indel_",VC, sep = "")) %>% 
  dplyr::rename(Tumour=Sample,gene=VD_Gene,variant_class=VC) %>% 
  left_join(id.key,by="Tumour") %>% relocate(sample,1) %>% 
  dplyr::select(c(sample,gene,variant_class)) %>% 
  filter(gene %in% top.TFs) %>% 
  distinct()

sub.drivers <- as.data.frame(read_excel("./data/SCANB/3_genomic/raw/HER2_enriched_June23_ForJohan.xlsx", sheet = "AllCodingASMD_CLPM")) %>% 
  mutate(VC = paste("sub_",VC, sep = "")) %>% 
  dplyr::rename(Tumour=Sample,gene=VD_Gene,variant_class=VC) %>% 
  left_join(id.key,by="Tumour") %>% relocate(sample,1) %>% 
  dplyr::select(c(sample,gene,variant_class)) %>% 
  filter(gene %in% top.TFs) %>% 
  distinct()

# processed amplification data
amp.dat <- loadRData("data/SCANB/4_CN/processed/CNA_genelevel.RData")[[2]] %>% 
  as.data.frame() %>% 
  filter(gene %in% top.TFs)

amp.dat <- pivot_longer(amp.dat,cols=colnames(amp.dat)[c(5:32)],
                        names_to = "sample",
                        values_to = "variant_class") %>% 
  dplyr::select(-c(chr,centerPos,Genome_pos)) %>% 
  filter(variant_class!=0) %>% mutate(variant_class="amplified")


mut.dat <- rbind(sub.drivers,amp.dat)

# make waterfall plot
plot <- waterfall(mut.dat, 
                  fileType = "Custom", 
                  variant_class_order = unique(mut.dat$variant_class),
                  mainGrid = TRUE,
                  main_geneLabSize = 15,
                  mainRecurCutoff = 0,
                  maxGenes = 30,
                  mainDropMut = FALSE, # drop unused mutation types from legend
                  #rmvSilent = TRUE,
                  out= "grob",
                  #mutBurdenLayer = layer,
                  plotMutBurden = FALSE)
#tbl <- read.csv("data/SCANB/3_genomic/raw/SCANBrel4_ExprSomMutations.csv")
grid::grid.newpage()
grid::grid.draw(plot)

################################################################################
# Calc. centroid over all 258 genes (or the two clusters) ERp-HER2E samples and correalte every sample to it (all ERp)-> make a heatmap of results to see if there are more samples belonging to that track
################################################################################

dev.off()