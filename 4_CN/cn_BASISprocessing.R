# Script: Process BASIS CN data to prepare for plotting (incl. Liftover)

# TODO:
# - convert count to % by dividing by # of samples
# - liftover of BASIS SNP-A positions to build 38 which SCAN-B is on.

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run (for plot naming)
cohort <- "BASIS" 

# set/create output directory for plots
output.path <- "output/plots/4_CN/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/4_CN/processed/",sep="")
dir.create(data.path)

#packages
source("scripts/4_CN/src/cn_functions.R")
library(ggplot2)
library(tidyverse)
library(readxl)
library(rtracklayer)
#BiocManager::install("liftOver")
library(liftOver)
library(GenomicRanges)
library(rtracklayer)
#library(R.utils)
library(Repitools)

################################################################################
# loading all CN data and convert GL freq. to %
################################################################################

# chr lengths
# get chr lengths to get genome positions of probes (excluding chr X)
chr.lengths <- as.data.frame(read.table(file = "data/BASIS/4_CN/raw/GRCh38_EBV.chrom.sizes.tsv", sep = '\t', header = FALSE))[1:22,] %>% 
  dplyr::rename(Chr=V1,length=V2) %>% 
  mutate(genome = cumsum(as.numeric(length))) %>% 
  mutate(genome = lag(genome,default = 0)) # lag by 1 position (cause I have to add the length of the previous chr to the probe positions (0 for chr1 probes))
chr.lengths$Chr <- as.numeric(gsub('^.{3}','',chr.lengths$Chr))

# calc GL alteration freqs in %
cn.luma <- loadRData("data/BASIS/4_CN/processed/LumA_CollectedFrequencyData.RData") 
cn.luma <- do.call("cbind", list(cn.luma$fData,cn.luma$CN_Gain,cn.luma$CN_Loss)) %>% 
  dplyr::rename(CN_Gain=ncol(.)-1,
                CN_Loss=ncol(.)) %>% 
  mutate(CN_Gain = (CN_Gain/length(cn.luma$Samples))*100) %>% 
  mutate(CN_Loss = (CN_Loss/length(cn.luma$Samples))*100)

cn.lumb <- loadRData("data/BASIS/4_CN/processed/Lumb_CollectedFrequencyData.RData") 
cn.lumb <- do.call("cbind", list(cn.lumb$fData,cn.lumb$CN_Gain,cn.lumb$CN_Loss)) %>% 
  dplyr::rename(CN_Gain=ncol(.)-1,
                CN_Loss=ncol(.)) %>% 
  mutate(CN_Gain = (CN_Gain/length(cn.lumb$Samples))*100) %>% 
  mutate(CN_Loss = (CN_Loss/length(cn.lumb$Samples))*100)

cn.her2p <- loadRData("data/BASIS/4_CN/processed/ERpHER2p_CollectedFrequencyData.RData") 
cn.her2p <- do.call("cbind", list(cn.her2p$fData,cn.her2p$CN_Gain,cn.her2p$CN_Loss)) %>% 
  dplyr::rename(CN_Gain=ncol(.)-1,
                CN_Loss=ncol(.)) %>% 
  mutate(CN_Gain = (CN_Gain/length(cn.her2p$Samples))*100) %>% 
  mutate(CN_Loss = (CN_Loss/length(cn.her2p$Samples))*100)

#View(head(cn.lumb))
#View(cn.her2p)
# lumb.samples <-cn.lumb$Samples
# lumb.samples
# length(lumb.samples)
# luma.samples <-cn.luma$Samples
# luma.samples
# length(luma.samples)
# 
# save(lumb.samples, 
#      file = "data/BASIS/3_genomic/processed/BASIS_LumB_samples.RData")
# save(luma.samples, 
#      file = "data/BASIS/3_genomic/processed/BASIS_LumA_samples.RData")


# l <- list(cn.luma,cn.lumb,cn.her2p)
# for (i in 1:3) {
#     str(l[i])
# }

################################################################################
# liftover of BASIS SNP-A positions to build 38 which SCAN-B is on.
################################################################################

# get the right chain file
#url <- "https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz"
#destfile <- "data/BASIS/4_CN/raw/hg19xToHg38.over.chain.gz"
#download.file(url, destfile)
#gunzip(destfile)
ch <- import.chain("data/BASIS/4_CN/raw/hg19ToHg38.over.chain")

# get df in right format to convert to granges object (only have to do that once for all files)
df <- cn.luma %>% 
  mutate(end = centerPosition) %>% 
  dplyr::rename(start = centerPosition) %>% 
  mutate(chromosome = paste("chr",as.character(chromosome),sep = "")) %>% 
  mutate(chromosome = if_else(chromosome == "chr23","chrX",chromosome)) %>% # change name to same as in chain file
  dplyr::select(-c(CN_Gain,CN_Loss))

# make granges object
gr.df <- makeGRangesFromDataFrame(df,seqnames.field = "chromosome",start.field = "start",end.field = "end", keep.extra.columns = TRUE,ignore.strand = TRUE) #starts.in.df.are.0based=FALSE # think not

# liftover
hg38.pos <- liftOver(gr.df,ch) 
hg38.pos <- unlist(hg38.pos) # convert to granges object
hg38.pos.df <- annoGR2DF(hg38.pos) %>% dplyr::select(-c("end","width")) 
hg38.pos.df$chr <- as.numeric(gsub("X",23,gsub('^.{3}','',hg38.pos.df$chr)))

# compare to see if the coordinates changed -> it worked :D
#annoGR2DF(gr.df) %>% filter(reporterId == "SNP_A-8575395") # 564477
#hg38.pos.df %>% filter(reporterId == "SNP_A-8575395") # 629097

################################################################################
# use liftover results to convert positions
################################################################################

cn.her2p.postlift <- merge(cn.her2p,hg38.pos.df,by="reporterId") %>% 
  dplyr::select(-c("chromosome","centerPosition")) %>% 
  dplyr::rename(Position = start, freqloss.HER2p = CN_Loss, 
                freqgain.HER2p = CN_Gain, Chr = chr, 
                ProbeID = reporterId) 

cn.luma.postlift <- merge(cn.luma,hg38.pos.df,by="reporterId") %>% 
  dplyr::select(-c("chromosome","centerPosition")) %>% 
  dplyr::rename(Position = start, freqloss.LUMA = CN_Loss, 
                freqgain.LUMA = CN_Gain, Chr = chr, 
                ProbeID = reporterId) 

cn.lumb.postlift <- merge(cn.lumb,hg38.pos.df,by="reporterId") %>% 
  dplyr::select(-c("chromosome","centerPosition")) %>% 
  dplyr::rename(Position = start, freqloss.LUMB = CN_Loss, 
                freqgain.LUMB = CN_Gain, Chr = chr, 
                ProbeID = reporterId) 

# merge
cn.basis.subtypes <- merge(cn.her2p.postlift,
                           merge(cn.luma.postlift, cn.lumb.postlift, 
                                 by = c("ProbeID","Chr","Position")),
                           by = c("ProbeID","Chr","Position"))

head(cn.basis.subtypes)

################################################################################
# Make key file for which probes map to which genes
################################################################################

probe.positions <- cn.basis.subtypes[1:3] %>% 
  mutate(Chr = paste("chr",as.character(Chr),sep = "")) %>% 
  filter(Chr != "chr23") # no X chr in scanb data

# create granges objects
probes <- GRanges(seqnames = probe.positions$Chr,
                  ranges = IRanges(probe.positions$Position),
                  ProbeID = probe.positions$ProbeID)
genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene) # meta = entrez ID

# convert to hgnc symbols
ENTREZID2SYMBOL <- select(org.Hs.eg.db, mcols(genes)$gene_id, c("ENTREZID", "SYMBOL"))
stopifnot(identical(ENTREZID2SYMBOL$ENTREZID, mcols(genes)$gene_id))
mcols(genes)$SYMBOL <- ENTREZID2SYMBOL$SYMBOL
gene.anno <- as.data.frame(genes)
probe.anno <- as.data.frame(probes)

# which 
overlap.res <- findOverlaps(genes,probes)
# queryHits(): indexes of the gene coordinates that overlap the corresponding 
# subjectHits(): indexes of the probes
# line up the query column identifier (gene) that overlaps each probe
f1 <- factor(subjectHits(overlap.res), levels=seq_len(subjectLength(overlap.res)))
# use of factor() with exactly as many levels as there are subjects ensures that the splitAsList() command returns a 1:1 mapping between the subjects (probes) and the genes in the corresponding CharacterList
overlap.list <- splitAsList(mcols(genes)[["SYMBOL"]][queryHits(overlap.res)], f1) # split the column of gene symbols into lists corresponding to the regions of overlap
mcols(probes) <- overlap.list

key.df <- merge(as.data.frame(probes),probe.anno,
                by=c("seqnames","start","end","width","strand")) %>% 
  dplyr::rename(Gene_symbol=X,Position=start,Chr=seqnames) %>% 
  dplyr::select(Chr,ProbeID,Gene_symbol)
#head(key.df) # why multiple genes per probe??

# add to matrix
probe.positions <- merge(probe.positions,key.df,by=c("Chr","ProbeID")) %>% 
  relocate(Gene_symbol,.after = ProbeID)

#head(probe.positions)

# replace character(0) with NA for later filtering of probes without annotation
probe.positions$Gene_symbol <- lapply(probe.positions$Gene_symbol, 
                                      function(x) {
                                        if (identical(x, character(0))) {
                                          return(NA)
                                        } else {return(x)}
                                      })

# filter NA rows
#probe.positions <- probe.positions[which(!is.na(probe.positions$Gene_symbol)),]

# which probes have multiple mapped genes
#multi.genes <- probe.positions[which(lengths(probe.positions$Gene_symbol)>1),]
#nrow(multi.genes) #107629

# key file for gene to probe mapping
map.key <- probe.positions

#View(map.key) 
  
################################################################################
# add genome position of probes and corresponding genes
################################################################################
# key file for gene to probe mapping - make one for BASIS or just use table14? 
#map.key <- loadRData("data//4_CN/processed/CN_mapped_probes.RData")

cnfile <- cn.basis.subtypes
# add genome position
cnfile <- add_genomepos(cnfile,chr.lengths)

# add "chr"
cnfile$Chr <- sub("^", "chr", cnfile$Chr) 

# merge with gene positions
cnfile <- merge(cnfile,map.key,by=c("Chr","ProbeID","Position")) %>% 
  relocate(Gene_symbol,.after = ProbeID)

# save
save(cnfile, 
     file = paste(data.path,"CN_gainloss_frequencies_genpos_genmap.RData",sep=""))


