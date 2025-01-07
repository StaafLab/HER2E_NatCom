# Script: create gene level CN matrices for SCANB and BASIS for later comparison
# based on the segment file with assigned CN status

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
library(purrr)
library(readxl)
library(IRanges)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(Repitools)
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
scanb.segments <- "data/SCANB/4_CN/processed/Segment_CN_states.RData"
basis.segments.1 <- "data/BASIS/4_CN/raw/ASCAT_CEL_Total/ASCAT_CEL_Total_EasySegments.RData"
basis.segments.2 <- "data/BASIS/4_CN/raw/ASCAT_InSilico/ASCAT_InSilico_EasySegments.RData"
chr.lengths <- "data/BASIS/4_CN/raw/GRCh38_EBV.chrom.sizes.tsv"
# output
plot.file <- "output/plots/4_CN/COMBINED_HER2n_signifgenes.pdf"
txt.file <- "output/plots/4_CN/COMBINED_HER2n_signifgenes.txt"
scanb.gene.cna <- "data/SCANB/4_CN/processed/CNA_genelevel.RData"
basis.gene.cna <- "data/BASIS/4_CN/processed/CNA_genelevel_all.RData" # del _all
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

# chr lengths
# get chr lengths to get genome positions of probes (excluding chr X)
chr.lengths <- as.data.frame(read.table(file = chr.lengths, sep = '\t', header = FALSE))[1:22,] %>% 
  dplyr::rename(chr=V1,length=V2) %>% 
  mutate(genome = cumsum(as.numeric(length))) %>% 
  mutate(genome = lag(genome,default = 0)) # lag by 1 position (cause I have to add the length of the previous chr to the probe positions (0 for chr1 probes))
chr.lengths$chr <- as.numeric(gsub('^.{3}','',chr.lengths$chr))

################################################################################
# SCANB get CN data on gene level (1 row per gene)
################################################################################

# get gene positions & convert to hgnc symbols
genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene) # meta = entrez ID)
ENTREZID2SYMBOL <- biomaRt::select(org.Hs.eg.db, mcols(genes)$gene_id, c("ENTREZID", "SYMBOL"))
stopifnot(identical(ENTREZID2SYMBOL$ENTREZID, mcols(genes)$gene_id))
mcols(genes)$SYMBOL <- ENTREZID2SYMBOL$SYMBOL
# convert to df and filter out NA gene symbols
genes <- annoGR2DF(genes) %>% 
  filter(!is.na(SYMBOL)) %>% 
  mutate(chr = gsub("^chr","",chr)) %>% 
  filter(chr %in% c(1:22)) %>% 
  dplyr::select(-c(gene_id,strand, width)) %>% 
  mutate_at(c('chr','start','end'), as.numeric)

# get segment data
cn.scanb.segments <- loadRData(scanb.segments)

# result storing objects
gl.names <- c("gene", "chr", "start", "end",
              c(names(cn.scanb.segments)))
amp.names <- c("gene", "chr", "start", "end",
               c(names(cn.scanb.segments)))
gl.df <- as.data.frame(matrix(nrow=length(genes$SYMBOL),ncol=length(gl.names)))
names(gl.df) <- gl.names
amp.df <- as.data.frame(matrix(nrow=length(genes$SYMBOL),ncol=length(amp.names)))
names(amp.df) <- amp.names

# loop over genes
pb = txtProgressBar(min = 0, max = length(genes$SYMBOL), initial = 0, style = 3)
for (i in 1:nrow(genes)) {
  setTxtProgressBar(pb,i)
  
  gene.dat <- genes[i,]

  # get the gainloss and amp state of that gene for all samples
  gene.statuses <- sapply(cn.scanb.segments, function(segment.df) {
    
    segment.df <- segment.df %>% 
      mutate_at(c('chr','startpos','endpos','GainLoss','Amp'), as.numeric)
    
    # get the cn status of the segment the gene is on
    # only relevant segments
    chr.segments <- segment.df[segment.df$chr==gene.dat$chr,]
    
    # check gene/segment overlap
    segments.query <- with(chr.segments, IRanges(startpos, endpos))
    gene.subject <- with(gene.dat, IRanges(start, end))
    
    # check 
    chr.segments$overlap = countOverlaps(segments.query, gene.subject) != 0
    if (sum(chr.segments$overlap)==0) { # no segment covers the gene region
      gene.gl.state <- NA
      gene.amp.state <- NA
    } else {
      hits <- findOverlaps(segments.query, gene.subject)
      overlaps <- pintersect(segments.query[queryHits(hits)], gene.subject[subjectHits(hits)])
      percentOverlap <- width(overlaps) / width(gene.subject[subjectHits(hits)])
      chr.segments <- chr.segments[chr.segments$overlap == TRUE,]
      chr.segments$percentOverlap <- percentOverlap
      # only select the segment with the highest overlap proportion
      main.segment <- chr.segments[which.max(chr.segments$percentOverlap),]
      gene.gl.state <- main.segment$GainLoss
      gene.amp.state <- main.segment$Amp
      }
    
    return(list("GainLoss"=gene.gl.state,"Amp"=gene.amp.state))
    })

  # store results

  gl.df[i,] <- c(gene.dat$SYMBOL,
                 gene.dat$chr,
                 gene.dat$start,
                 gene.dat$end,
                 unname(unlist(gene.statuses["GainLoss",])))
  
  amp.df[i,] <- c(gene.dat$SYMBOL,
                  gene.dat$chr,
                  gene.dat$start,
                  gene.dat$end,
                  unname(unlist(gene.statuses["Amp",])))

  close(pb)
}

cn.list <- list("gainloss"=gl.df,"amp"=amp.df)

# add center position
cn.list <- lapply(cn.list, function(df) {
  df <- df %>% mutate_at(c("chr","start","end"), as.numeric)
  centerPos.vec <- apply(df, 1, function(row) {
    gene.length <- as.numeric(row[["end"]]) - as.numeric(row[["start"]])
    centerPos <- as.numeric(row[["start"]]) + (gene.length/2)
    return(centerPos)
  })
  df$centerPos <- centerPos.vec
  df <- df %>% dplyr::relocate(centerPos, .after=chr)
  return(df)
})

# convert to genome position
cn.list <- lapply(cn.list, function(df) {
  df <- df %>% 
    mutate_at(c("chr","start","end"), as.numeric) %>% 
    # add new chr genome position column
    mutate(genome = 0) %>% 
    # update the genome col to fill in the actual chr positions
    rows_update(chr.lengths[c("chr","genome")]) %>% 
    # add a column with the genome position of each probe
    mutate(Genome_pos = centerPos + genome) %>% 
    relocate(c(genome,Genome_pos), .after=centerPos) %>% 
    dplyr::select(-c(genome))
  return(df)
})

# remove start and end
cn.list <- lapply(cn.list, function(df) {
  df <- df %>% 
    dplyr::select(-c(start, end))
  return(df)
})

#save(cn.list,
#     file = scanb.gene.cna)

cn.scanb.list <- loadRData(scanb.gene.cna)

################################################################################
# BASIS get CN data on gene level (1 row per gene)
################################################################################

# hg19 genes
# get gene positions & convert to hgnc symbols
genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
ENTREZID2SYMBOL <- biomaRt::select(org.Hs.eg.db, mcols(genes)$gene_id, c("ENTREZID", "SYMBOL"))
stopifnot(identical(ENTREZID2SYMBOL$ENTREZID, mcols(genes)$gene_id))
mcols(genes)$SYMBOL <- ENTREZID2SYMBOL$SYMBOL
# filter only genes found in the scanb file
genes <- annoGR2DF(genes) %>% 
  filter(!is.na(SYMBOL)) %>% 
  mutate(chr = gsub("^chr","",chr)) %>% 
  filter(chr %in% c(1:22)) %>% 
  dplyr::select(-c(gene_id,strand, width)) %>% 
  mutate_at(c('chr','start','end'), as.numeric) %>% 
  filter(SYMBOL %in% cn.scanb.list[[1]]$gene)

# get segment data
cn.basis.segments.1 <- loadRData(basis.segments.1)
cn.basis.segments.2 <- loadRData(basis.segments.2)
cn.basis.segments <- rbind(cn.basis.segments.1,cn.basis.segments.2)

# get relevant sample IDs
basis.samples <- as.data.frame(loadRData("data/BASIS/1_clinical/raw/Summarized_Annotations_BASIS.RData")) %>% 
  dplyr::select(c(sample_name,final.ER,final.HER2,
                  ClinicalGroup,PAM50_AIMS,MutSigProportions)) %>%
  mutate(Subtype = 
           case_when(final.ER=="positive" & final.HER2=="negative" & PAM50_AIMS=="LumA" ~ "LumA",
                     final.ER=="positive" & final.HER2=="negative" & PAM50_AIMS=="LumB" ~ "LumB")) %>% 
  dplyr::filter(Subtype %in% c("LumA","LumB")) %>% 
  pull(sample_name)

# select relevant samples
cn.basis.segments$SampleID <- gsub("b","",gsub("a.*","",cn.basis.segments$SampleID))
cn.basis.segments <- cn.basis.segments %>% filter(SampleID %in% basis.samples)

# organize in list of dfs (1 per sample)
cn.basis.segments <- with(cn.basis.segments, split(
  cn.basis.segments, list(SampleID)))
cn.basis.segments <- lapply(cn.basis.segments, function(x) {
  x <- x %>% dplyr::rename(GainLoss=CNA,startpos=start,endpos=stop)
  return(x)
})

# no amp data for basis

# result storing objects
gl.names <- c("gene", "chr", "start", "end",
              c(names(cn.basis.segments)))
gl.df <- as.data.frame(matrix(nrow=length(genes$SYMBOL),ncol=length(gl.names)))
names(gl.df) <- gl.names

# loop over genes
start.time <- Sys.time()
pb = txtProgressBar(min = 0, max = length(genes$SYMBOL), initial = 0, style = 3)
for (i in 1:nrow(genes)) {
  setTxtProgressBar(pb,i)
  
  gene.dat <- genes[i,]
  
  # get the gainloss state of that gene for all samples
  gene.statuses <- sapply(cn.basis.segments, function(segment.df) {
    segment.df <- segment.df %>% 
      mutate_at(c('chr','startpos','endpos','GainLoss'), as.numeric)
    # get the cn status of the segment the gene is on
    # only relevant segments
    chr.segments <- segment.df[segment.df$chr==gene.dat$chr,]
    
    segments.query <- with(chr.segments, IRanges(startpos, endpos))
    gene.subject <- with(gene.dat, IRanges(start, end))
    # check 
    chr.segments$overlap = countOverlaps(segments.query, gene.subject) != 0
    if (sum(chr.segments$overlap)==0) { # no segment covers the gene region
      gene.gl.state <- NA
    } else {
      hits <- findOverlaps(segments.query, gene.subject)
      overlaps <- pintersect(segments.query[queryHits(hits)], gene.subject[subjectHits(hits)])
      percentOverlap <- width(overlaps) / width(gene.subject[subjectHits(hits)])
      chr.segments <- chr.segments[chr.segments$overlap == TRUE,]
      chr.segments$percentOverlap <- percentOverlap
      # only select the segment with the highest overlap proportion
      main.segment <- chr.segments[which.max(chr.segments$percentOverlap),]
      gene.gl.state <- main.segment$GainLoss
    }
    
    return(list("GainLoss"=gene.gl.state,"Amp"=NA))
  })
  
  # store results
  gl.df[i,] <- c(gene.dat$SYMBOL,
                 gene.dat$chr,
                 gene.dat$start,
                 gene.dat$end,
                 unname(unlist(gene.statuses["GainLoss",])))
  
  close(pb)
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# dont need to do add center position, convert to genome position, and liftover as I can just take that data from scanb genes
scanb.gl <- cn.scanb.list[[1]]
gl.df$centerPos <- scanb.gl$centerPos[match(gl.df$gene,scanb.gl$gene)]
gl.df$Genome_pos <- scanb.gl$Genome_pos[match(gl.df$gene,scanb.gl$gene)]
gl.df$chr <- scanb.gl$chr[match(gl.df$gene,scanb.gl$gene)]
gl.df <- gl.df %>% 
  dplyr::select(-c(start,end)) %>% 
  relocate(c(chr,centerPos,Genome_pos), .after=gene)

#save(gl.df, 
#     file = basis.gene.cna)
gl.df <- loadRData(basis.gene.cna)

View(gl.df)
