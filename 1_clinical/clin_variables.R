# Script: Create table with clinicopathological variables comparisons between subtypes

# if run for SCANB, output: 
# - table whole rel 4 + h2low +% review
# - table review 1 

# for METABRIC:
# - table + IC10

# note: I store the sd for continuous vars in the % column

# TODO: - rerun with new RS1 sample set (trement filtered) and add SD to age and size!

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "Metabric" # Metabric or SCANB 

# set/create output directory for plots
output.path <- "output/plots/1_clinical/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/1_clinical/processed/",sep="")
dir.create(data.path)

#packages
source("scripts/1_clinical/src/clin_functions.R")
library(tidyverse)
library(matrixStats)
library(DescTools)
library(openxlsx)
library(writexl)

source_dat_path <- "./output/source_data/R_objects/"  

#######################################################################
# Cohort-specific data preprocessing including selection of  
# the clinical ER+Her2- subtyped samples
#######################################################################

# for metabric cohort
if (cohort=="Metabric") {
    
    # load data
    # extract relevant variables
    anno <- loadRData("./data/METABRIC/1_clinical/processed/Merged_annotations_ERpHER2n.RData") %>% 
        mutate(HER2_Low = case_when(
          HER2_IHC_status==1 ~ 1,
          HER2_IHC_status==2 & HER2_SNP6_state!="GAIN" ~ 1)) %>% 
        mutate(LN = ifelse(lymph_nodes_positive > 0, 1, 0)) %>%  
        dplyr::rename(sampleID=METABRIC_ID,
                      NHG=Grade,
                      Size=TumSize,
                      IC10=IC10_Rueda) %>% 
        dplyr::select(sampleID, PAM50, Age, Size, NHG, LN, HER2_Low, IC10) # no pr status available
    
    # replace na values with 0
    anno$HER2_Low[is.na(anno$HER2_Low)] <- 0
    
    # 
    anno$IC10 <- factor(anno$IC10, levels=c("1","2","3","4ER+","6","7","8","9","10"))
    # save source data
    save(anno, file=paste0(source_dat_path,"Table_1_Metabric.RData")) # Save 
    
# for SCANB cohort
} else if (cohort=="SCANB") {
    
    # select relevant variables
    anno <- loadRData("./data/SCANB/1_clinical/processed/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data_ERpHER2n.RData") %>% 
        mutate(LN_Bosch = case_when(LNstatus_Bosch == "N0" ~ 0, LNstatus_Bosch == "N+" ~ 1)) %>% 
      dplyr::select(
        GEX.assay, 
        # table 1: N, Age, size, PR, NHG, LN, % review, H2low
        NCN.PAM50, Age, PR, Size.mm, NHG, LN, HER2_Low,
        # table 2: N, Age, size, PR, NHG, LN, H2low
        Bosch_RS1, TumSize_Bosch, NHG_Bosch, LN_Bosch
        ) %>% 
      dplyr::rename(
            sampleID=GEX.assay,
            PAM50=NCN.PAM50,
            Size=Size.mm,
            Size_Bosch = TumSize_Bosch)
    
    # do data type correction for scanb exclusive vars 
    #anno$PR <- as.factor(anno$PR)
    anno$Size_Bosch <- as.numeric(anno$Size_Bosch)
    anno$NHG_Bosch <- as.factor(anno$NHG_Bosch)
    anno$LN_Bosch <- as.factor(anno$LN_Bosch)
    anno$PR[anno$PR == ""] <- NA
    
    # save source data
    anno.export <- anno
    names(anno.export)[names(anno.export) == "Bosch_RS1"] <- "Review_cohort"
    colnames(anno.export) <- gsub("_Bosch", "_review", colnames(anno.export))
    save(anno.export, file=paste0(source_dat_path,"Table_1_SCANB.RData")) # Save 
    
}

# get correct data types
anno$PAM50 <- as.factor(anno$PAM50)
anno$NHG[anno$NHG == ""] <- NA
anno$NHG <- as.factor(anno$NHG)
anno$Size <- as.numeric(anno$Size)
anno$Age <- as.numeric(anno$Age)
anno$LN <- as.factor(anno$LN)
anno$HER2_Low <- as.factor(anno$HER2_Low)

# check if i can do that differently
# sample names
Her2.ids <- anno %>% filter(PAM50=="Her2") %>% pull(sampleID)
LumA.ids <- anno %>% filter(PAM50=="LumA") %>% pull(sampleID)
LumB.ids <- anno %>% filter(PAM50=="LumB") %>% pull(sampleID)
anno <- anno %>% column_to_rownames(var = "sampleID")

# storing format: excel document with sheets for each variable
# keep column names consitent so later I just need to copy & paste tables together
col.names <- c("Variable","HER2E(ref)","HER2E.%","LUMA","LUMA.%","LUMA.pval","LUMB","LUMB.%","LUMB.pval")

########################################################################
# 1. data vars (not review)
# table 1: N, Age, size, PR, NHG, LN, % review, H2low

#######################################################################
# PAM50 count
#######################################################################
# matrix
pam50.count <- as.data.frame(matrix(ncol=length(col.names)))
names(pam50.count) <- col.names
pam50.count$Variable <- "N"

t <- table(anno$PAM50)

# add count data
pam50.count$`HER2E(ref)`[1] <- t[1]
pam50.count$LUMA[1] <- t[2]
pam50.count$LUMB[1] <- t[3]

# add % data
pam50.count$`HER2E.%`[1] <- (t[1]/length(anno$PAM50))*100
pam50.count$`LUMA.%`[1] <- (t[2]/length(anno$PAM50))*100
pam50.count$`LUMB.%`[1] <- (t[3]/length(anno$PAM50))*100

pam50.count

#######################################################################
# Age comparison
#######################################################################
# note: I store the sd for continuous vars in the % column

age.data <- anno %>% filter(!is.na(Age))
table(anno$PAM50,is.na(anno$Age))

# matrix
age.df <- as.data.frame(matrix(ncol=length(col.names)))
names(age.df) <- col.names
age.df$Variable <- "Age"

age.df <- test.cont(data = age.data,
                    var.df = age.df,
                    var = "Age")

age.df

#######################################################################
# 4. Size
#######################################################################

size.data <- anno %>% filter(!is.na(Size))
table(anno$PAM50,is.na(anno$Size))

# matrix
size.df <- as.data.frame(matrix(ncol=length(col.names)))
names(size.df) <- col.names
size.df$Variable <- "Size"

size.df <- test.cont(data = size.data,
          var.df = size.df,
          var = "Size")

size.df

#######################################################################
# 4. NHG (2x3 contingency table)
#######################################################################

nhg.data <- anno %>% filter(!is.na(NHG))
table(anno$PAM50,is.na(anno$NHG))

# matrix
nhg.df <- as.data.frame(matrix(ncol=length(col.names),nrow = 3))
names(nhg.df) <- col.names
nhg.df$Variable <- c("NHG1","NHG2","NHG3")

nhg.df <-   test.2x3ct(data = nhg.data,
                       var.df = nhg.df,
                       var = "NHG")

nhg.df

#######################################################################
# 4. LN (2x2 contingency table)
#######################################################################

ln.data <- anno %>% filter(!is.na(LN))
table(anno$PAM50,is.na(anno$LN))

# matrix
ln.df <- as.data.frame(matrix(ncol=length(col.names),nrow = 2))
names(ln.df) <- col.names
ln.df$Variable <- c("LN.N0","LN.N+")

ln.df <- test.2x2ct(data = ln.data,
           var.df = ln.df,
           var = "LN")

#######################################################################
# 4. IC10 (only Metabric)
if (cohort=="Metabric") {
#######################################################################
  
  var <- "IC10"
  data <- anno %>% filter(!is.na(IC10))
  table(anno$PAM50,is.na(anno$IC10))
  
  # matrix
  var.df <- as.data.frame(matrix(ncol=length(col.names),nrow = 9))
  names(var.df) <- col.names
  var.df$Variable <- c("IC-1","IC-2","IC-3","IC-4ER+","IC-6","IC-7","IC-8","IC-9","IC-10") # no ic5 classification, add later in table
  
  # for Her2 vs LumA
  luma.ct <- table(data$PAM50,data[[var]])[c("Her2","LumA"),]
    
  # for Her2 vs LumB
  lumb.ct <- table(data$PAM50,data[[var]])[c("Her2","LumB"),]
  
  groups <- c("Her2","LumA","LumB")
  cols <- c("HER2E","LUMA","LUMB")
  
  for (i in 1:3) {
    type <- groups[i]
    col <- cols[i]
    type.dat <- data[which(data$PAM50==type),]
    type.dat$PAM50 <- droplevels(type.dat$PAM50)
    type.counts <- table(type.dat$PAM50,type.dat[[var]])
    
    # count column
    for (j in 1:length(unique(data$IC10))) {
      if (col=="HER2E") {
        var.df[[paste(col,"(ref)",sep="")]][j] <- type.counts[j]
      } else {
        var.df[[col]][j] <- type.counts[j]
      }
    
      # % column
      var.df[[paste(col,".%",sep="")]][j] <- round(type.counts[j]/sum(type.counts)*100)
      var.df[[paste(col,".%",sep="")]][j] <- round(type.counts[j]/sum(type.counts)*100)
      var.df[[paste(col,".%",sep="")]][j] <- round(type.counts[j]/sum(type.counts)*100)
    }
  }
  
  ic10.df <- var.df
  
  # create grouped barplot
  dat <- as.data.frame(table(anno$PAM50,anno$IC10)) %>% 
    mutate(Freq=case_when(Var1=="Her2" ~ (Freq/table(anno$PAM50)[["Her2"]])*100,
                          Var1=="LumA" ~ (Freq/table(anno$PAM50)[["LumA"]])*100,
                          Var1=="LumB" ~ (Freq/table(anno$PAM50)[["LumB"]])*100))
  
  # Create grouped barplot using ggplot2
  ic10.plot <- ggplot(dat,aes(x = Var2, y =Freq, fill = Var1)) +
    geom_bar(stat = "identity", color="black",width=0.7, position = "dodge") +
    scale_fill_manual(name = "PAM50 Subtype",values=c("#d334eb", "#2176d5","#34c6eb")) +
    xlab("IC10 Subtype") +
    ylab("Fraction (%)")
  
  #######################################################################
  # HER2low frequency (2x2 contingency table)
  #######################################################################
  
  h2low.data <- anno %>% filter(!is.na(HER2_Low))
  table(anno$PAM50,anno$HER2_Low)
  
  # matrix
  h2low.df <- as.data.frame(matrix(ncol=length(col.names),nrow = 2))
  names(h2low.df) <- col.names
  h2low.df$Variable <- c("HER2low.No","HER2low.Yes")
  
  h2low.df <- test.2x2ct(data = h2low.data,
                         var.df = h2low.df,
                         var = "HER2_Low")
  
  
  h2low.df
}

#######################################################################
# only for SCANB
if (cohort=="SCANB") {
  #######################################################################
  # 4. PR status (2x2 contingency table)
  #######################################################################
  
  pr.data <- anno %>% filter(!is.na(PR))
  
  # matrix
  pr.df <- as.data.frame(matrix(ncol=length(col.names),nrow = 2))
  names(pr.df) <- col.names
  pr.df$Variable <- c("PR.neg","PR.pos")
  
  pr.df <- test.2x2ct(data = pr.data,
                      var.df = pr.df,
                      var = "PR")
  pr.df
  
  
  ########################################################################
  # 2. review data vars
  
  #######################################################################
  # % reviewed 
  #######################################################################
  
  # matrix
  rev.count <- as.data.frame(matrix(ncol=length(col.names)))
  names(rev.count) <- col.names
  rev.count$Variable <- "N.Reviewed"
  
  t.rev <- table(anno[which(anno$Bosch_RS1==1),]$PAM50)
  t.all <- table(anno$PAM50)
  
  # add count data
  rev.count$`HER2E(ref)`[1] <- t.rev[1]
  rev.count$LUMA[1] <- t.rev[2]
  rev.count$LUMB[1] <- t.rev[3]
  
  # add % data
  rev.count$`HER2E.%`[1] <- (t.rev[1]/t.all[1])*100
  rev.count$`LUMA.%`[1] <- (t.rev[2]/t.all[2])*100
  rev.count$`LUMB.%`[1] <- (t.rev[3]/t.all[3])*100
  
  rev.count
  
  #######################################################################
  # filter anno to only include review data
  anno <- anno %>% filter(Bosch_RS1 == 1)
  
  #######################################################################
  # HER2low frequency (2x2 contingency table)
  #######################################################################
  
  h2low.data <- anno %>% filter(!is.na(HER2_Low))
  table(anno$PAM50,anno$HER2_Low)
  
  # matrix
  h2low.df <- as.data.frame(matrix(ncol=length(col.names),nrow = 2))
  names(h2low.df) <- col.names
  h2low.df$Variable <- c("HER2low.No","HER2low.Yes")
  
  h2low.df <- test.2x2ct(data = h2low.data,
                         var.df = h2low.df,
                         var = "HER2_Low")
  
  
  h2low.df
  
  #######################################################################
  # PAM50 count in review
  #######################################################################
  
  # matrix
  pam50.count.rev <- as.data.frame(matrix(ncol=length(col.names)))
  names(pam50.count.rev) <- col.names
  pam50.count.rev$Variable <- "N.Bosch"
  
  t <- table(anno$PAM50)
  
  # add count data
  pam50.count.rev$`HER2E(ref)`[1] <- t[1]
  pam50.count.rev$LUMA[1] <- t[2]
  pam50.count.rev$LUMB[1] <- t[3]
  
  # add % data
  pam50.count.rev$`HER2E.%`[1] <- (t[1]/length(anno$PAM50))*100
  pam50.count.rev$`LUMA.%`[1] <- (t[2]/length(anno$PAM50))*100
  pam50.count.rev$`LUMB.%`[1] <- (t[3]/length(anno$PAM50))*100
  
  pam50.count.rev
  
  #######################################################################
  # LN Bosch
  #######################################################################
  
  ln_bosch.data <- anno %>% filter(!is.na(LN_Bosch))
  table(anno$PAM50,anno$LN_Bosch)
  
  # matrix
  ln_bosch.df <- as.data.frame(matrix(ncol=length(col.names),nrow = 2))
  names(ln_bosch.df) <- col.names
  ln_bosch.df$Variable <- c("LN.N0.Bosch","LN.N+.Bosch")
  
  ln_bosch.df <- test.2x2ct(data = ln_bosch.data,
                            var.df = ln_bosch.df,
                            var = "LN_Bosch")
  
  ln_bosch.df
  
  #######################################################################
  # Size Bosch x.data .df VAR
  #######################################################################
  
  size_bosch.data <- anno %>% filter(!is.na(Size_Bosch))
  table(anno$PAM50,is.na(anno$Size_Bosch))
  
  # matrix
  size_bosch.df <- as.data.frame(matrix(ncol=length(col.names)))
  names(size_bosch.df) <- col.names
  size_bosch.df$Variable <- "Size.Bosch"
  
  size_bosch.df <- test.cont(data = size_bosch.data,
                             var.df = size_bosch.df,
                             var = "Size_Bosch")
  
  size_bosch.df
  
  #######################################################################
  # NHG Bosch
  #######################################################################
  
  nhg_bosch.data <- anno %>% filter(!is.na(NHG_Bosch))
  table(anno$PAM50,is.na(anno$NHG_Bosch))
  
  # matrix
  nhg_bosch.df <- as.data.frame(matrix(ncol=length(col.names),nrow = 3))
  names(nhg_bosch.df) <- col.names
  nhg_bosch.df$Variable <- c("NHG1.Bosch","NHG2.Bosch","NHG3.Bosch")
  
  nhg_bosch.df <-   test.2x3ct(data = nhg_bosch.data,
                         var.df = nhg_bosch.df,
                         var = "NHG_Bosch")
  
  nhg_bosch.df
  
  #######################################################################
  # Age for review
  #######################################################################
  
  age_bosch.data <- anno %>% filter(!is.na(Age))
  table(anno$PAM50,is.na(anno$Age))
  
  # matrix
  age_bosch.df <- as.data.frame(matrix(ncol=length(col.names)))
  names(age_bosch.df) <- col.names
  age_bosch.df$Variable <- "Age"
  
  age_bosch.df <- test.cont(data = age_bosch.data,
                      var.df = age_bosch.df,
                      var = "Age")
  
  age_bosch.df
  
  #######################################################################
  # PR for review
  #######################################################################
  
  pr_bosch.data <- anno %>% filter(!is.na(PR))
  
  # matrix
  pr_bosch.df <- as.data.frame(matrix(ncol=length(col.names),nrow = 2))
  names(pr_bosch.df) <- col.names
  pr_bosch.df$Variable <- c("PR.neg","PR.pos")
  
  pr_bosch.df <- test.2x2ct(data = pr_bosch.data,
                      var.df = pr_bosch.df,
                      var = "PR")
  pr_bosch.df
  
  #######################################################################
  # export to excel
  #######################################################################
    
  # rel4+part of rs1
  sheet.list.table1 <- list("N" = pam50.count, 
                           "Age" = age.df,
                           "PR" = pr.df,
                           "Size" = size.df,
                           "NHG" = nhg.df,
                           "LN" = ln.df,
                           "%reviewed" = rev.count,
                           "HER2low" = h2low.df)
  
  # rbind
  table1 <- as.data.frame(do.call("rbind", sheet.list.table1)) %>% mutate_at(vars(colnames(.)[colnames(.) %!in% c("Variable","LUMA.pval","LUMB.pval")]), function(x) {round(x,1)})
  
  write_xlsx(table1,"./output/supplementary_data/SCANB_clinpath_table1.xlsx")
  
  # review
  sheet.list.table2 <- list("N" = pam50.count.rev, 
                            "Age" = age_bosch.df,
                            "PR" = pr_bosch.df,
                            "Size" = size_bosch.df,
                            "NHG" = nhg_bosch.df,
                            "LN" = ln_bosch.df,
                            "HER2low" = h2low.df)
  
  # rbind
  table2 <- as.data.frame(do.call("rbind", sheet.list.table2)) %>% mutate_at(vars(colnames(.)[colnames(.) %!in% c("Variable","LUMA.pval","LUMB.pval")]), function(x) {round(x,1)})
  
  write_xlsx(table2,"./output/supplementary_data/SCANB_clinpath_review_table2.xlsx")
  
} else if (cohort=="Metabric") {
  # 1 table
  sheet.list.table1 <- list("N" = pam50.count, 
                            "Age" = age.df,
                            "Size" = size.df,
                            "NHG" = nhg.df,
                            "LN" = ln.df,
                            "HER2low" = h2low.df,
                            "IC10" = ic10.df) 
  
  # rbind
  table1 <- as.data.frame(do.call("rbind", sheet.list.table1)) %>% mutate_at(vars(colnames(.)[colnames(.) %!in% c("Variable","LUMA.pval","LUMB.pval")]), function(x) {round(x,1)})
  
  write_xlsx(table1,"./output/supplementary_data/Metabric_clinpath_table1.xlsx")
  
  # save plot
  pdf(file = paste(output.path,"Metabric_IC10_barplot.pdf",sep=""), 
      width = 21/2, height = 14.8/2) 
  print(ic10.plot)
  dev.off()
}

