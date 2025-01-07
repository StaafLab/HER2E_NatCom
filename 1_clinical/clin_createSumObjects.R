# Script: create summary objects for SCANB and Metabric
#TODO: - change the BOSCH rs1 column to 0 for cases that doint have review treatment annotation

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

#packages
source("./scripts/1_clinical/src/clin_functions.R")
library(ggplot2)
library(tidyverse)
library(matrixStats)
library(pheatmap)
library(readxl)
library(ggfortify)
library(janitor)
library(biomaRt)

##########################################################################
# create summary objects
# for SCANB review: exclude cases that have no treatment anno
##########################################################################

# 1. SCANB ; ERpHER2n ; LUMA ; LUMB ; HER2E
SCANB.ERpHER2n.data <- loadRData(
  "./data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData") %>% 
  filter(Follow.up.cohort==TRUE) %>% 
  filter(fuV8==TRUE) %>% 
  filter(ER=="Positive") %>% 
  filter(HER2=="Negative") %>% 
  filter(NCN.PAM50 %in% c("LumA","LumB","Her2")) %>% 
  mutate(TreatmentGroup_Bosch = case_when(
    Chemo_Bosch & ET_Bosch ~ "CE", 
    !Chemo_Bosch & ET_Bosch ~ "E")) %>% 
  mutate(Bosch_RS1 = ifelse(!is.na(TreatmentGroup_Bosch), 1, NA))

# add the metagene scores
SCANB.mg.scores <- loadRData(
  "./data/SCANB/2_transcriptomic/processed/mg_anno.RData")[[1]] %>% 
  dplyr::select(-c(PAM50)) %>% 
  dplyr::rename(GEX.assay = "sampleID")
SCANB.ERpHER2n.mg.data <- merge(SCANB.ERpHER2n.data,SCANB.mg.scores,by="GEX.assay") #,all=TRUE

#table(SCANB.ERpHER2n.mg.data$NCN.PAM50,is.na(SCANB.ERpHER2n.mg.data$Basal))

save(SCANB.ERpHER2n.mg.data,file="./data/SCANB/1_clinical/processed/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data_ERpHER2n.RData")

#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#

# 2. SCANB ; ERpHER2p ; ERpHER2nHER2E
SCANB.ERpHER2p.data <- loadRData(
  "./data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData") %>% 
  filter(Follow.up.cohort==TRUE) %>% 
  filter(fuV8==TRUE) %>% 
  filter(ER=="Positive") %>% 
  mutate(Group = case_when(
    HER2 == "Negative" & NCN.PAM50 == "Her2" ~ "HER2n_HER2E",
    HER2 == "Positive" & NCN.PAM50 == "Her2" ~ "HER2p_HER2E",
    HER2 == "Positive" & NCN.PAM50 != "Her2" ~ "HER2p_nonHER2E")) %>% 
  filter(Group %in% 
           c("HER2n_HER2E","HER2p_HER2E","HER2p_nonHER2E"))

# add metagene scores
SCANB.mg.scores <- loadRData(
  "./data/SCANB/2_transcriptomic/processed/mg_anno_HER2p.RData")[[1]] %>% 
  dplyr::select(-c(Group)) %>% 
  dplyr::rename(GEX.assay = "sampleID")
SCANB.ERpHER2p.mg.data <- merge(SCANB.ERpHER2p.data,SCANB.mg.scores,by="GEX.assay") #,all=TRUE

#table(SCANB.ERpHER2p.mg.data$NCN.PAM50,is.na(SCANB.ERpHER2p.mg.data$Basal))

save(SCANB.ERpHER2p.mg.data,file="./data/SCANB/1_clinical/processed/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data_ERpHER2p.RData")

#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#

# 3. METARBRIC ; ERpHER2n ; LUMA ; LUMB ; HER2E
#just cheking definition of ERpHER2n
#gg <- anno[which(anno$ER_IHC_status=="pos" & anno$HER2_SNP6_state!="GAIN" & anno$HER2_SNP6_state!="UNDEF"),]
#setequal(gg$METABRIC_ID,anno[which(anno$ClinGroup == "ERpHER2n"),]$METABRIC_ID)

METABRIC.ERpHER2n.data <- loadRData(
  "data/METABRIC/1_clinical/raw/Merged_annotations.RData") %>% 
  filter(PAM50 %in% c("LumA", "LumB", "Her2")) %>% 
  filter(grepl('ERpHER2n', ClinGroup))

METABRIC.mg.scores <- loadRData(
  "./data/METABRIC/2_transcriptomic/processed/mg_anno.RData")[[1]] %>% 
  dplyr::select(-c(PAM50)) %>% 
  dplyr::rename(METABRIC_ID = "sampleID")
METABRIC.ERpHER2n.mg.data <- merge(METABRIC.ERpHER2n.data,METABRIC.mg.scores,by="METABRIC_ID") #,all=TRUE

save(METABRIC.ERpHER2n.mg.data,file="./data/METABRIC/1_clinical/processed/Merged_annotations_ERpHER2n.RData")

#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#

# 4. METABRIC ; ERpHER2p ; ERpHER2nHER2E
# extract relevant variables
METABRIC.ERpHER2p.data <- loadRData(
  "data/METABRIC/1_clinical/raw/Merged_annotations.RData") %>% 
  mutate(ER = if_else(
    ER_IHC_status=="pos","Positive","Negative")) %>% 
  mutate(HER2 = case_when(
    HER2_SNP6_state == "GAIN" ~ "Positive",
    HER2_SNP6_state == "UNDEF" ~ "Undefined",
    HER2_SNP6_state %in% c("LOSS","NEUT") ~ "Negative")) %>% 
  filter(ER == "Positive") %>% 
  dplyr::select(-c("group")) %>% 
  mutate(Group = case_when(
    HER2 == "Negative" & PAM50 == "Her2" ~ "HER2n_HER2E",
    HER2 == "Positive" & PAM50 == "Her2" ~ "HER2p_HER2E",
    HER2 == "Positive" & PAM50 != "Her2" ~ "HER2p_nonHER2E")) %>% 
  filter(Group %in% 
           c("HER2n_HER2E","HER2p_HER2E","HER2p_nonHER2E"))

METABRIC.mg.scores <- loadRData(
  "./data/METABRIC/2_transcriptomic/processed/mg_anno_HER2p.RData")[[1]] %>% 
  dplyr::select(-c(Group)) %>% 
  dplyr::rename(METABRIC_ID = "sampleID")
METABRIC.ERpHER2p.mg.data <- merge(METABRIC.ERpHER2p.data,METABRIC.mg.scores,by="METABRIC_ID")

#table(METABRIC.ERpHER2p.mg.data$PAM50,is.na(METABRIC.ERpHER2p.mg.data$Basal))

save(METABRIC.ERpHER2p.mg.data,file="./data/METABRIC/1_clinical/processed/Merged_annotations_ERpHER2p.RData")


#
# load("./data/SCANB/1_clinical/processed/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data_ERpHER2n.RData")
# load("./data/SCANB/1_clinical/processed/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data_ERpHER2p.RData")
# load("./data/METABRIC/1_clinical/processed/Merged_annotations_ERpHER2n.RData")
# load("./data/METABRIC/1_clinical/processed/Merged_annotations_ERpHER2p.RData")
# 
# colnames(SCANB.ERpHER2n.mg.data)
# colnames(SCANB.ERpHER2p.mg.data)
# colnames(METABRIC.ERpHER2n.mg.data)
# colnames(METABRIC.ERpHER2p.mg.data)
# table(SCANB.ERpHER2n.mg.data$NCN.PAM50)
# table(SCANB.ERpHER2p.mg.data$Group)
# table(METABRIC.ERpHER2n.mg.data$PAM50)
# table(METABRIC.ERpHER2p.mg.data$Group)
