# Function definitions for transcriptomic data analyses

library(rstatix)
library(stats)

################################################################################
# functions
################################################################################

# loads RData data file and allows to assign it directly to variable
loadRData <- function(file.path){
    load(file.path)
    get(ls()[ls() != "file.path"])
}

# preprocessing scanb gex data to have hgnc as the gene ids
scanb_gex_load <- function(gex.path,geneanno.path,ID.type) {
    
    # load gex data
    gex.data <- as.data.frame(loadRData(gex.path))
    
    # load gene anno data to convert IDs
    gene.anno <- as.data.frame(loadRData(geneanno.path))
    
    # process
    gex.data <- gex.data %>% 
        rownames_to_column("ensembl_gene_id") %>% 
        mutate(geneID = gene.anno[which(gene.anno$Gene.ID==ensembl_gene_id),][[ID.type]]) %>% #gene.name=Idtype
        dplyr::select(-c(ensembl_gene_id)) %>% 
        drop_na(geneID) %>% 
        distinct(geneID,.keep_all = TRUE) %>% 
        column_to_rownames("geneID")
    
    return(gex.data)
    
}

# preprocessing metabric gex data to have hgnc as the gene ids
metabric_gex_load <- function(gex.path,ID.type) {
    
    # load gex data
    gex.data <- as.data.frame(read.table(gex.path, sep="\t"))
    
    # the other gene id to remove
    if(ID.type=="Hugo_Symbol") {
        ID.remove <- "Entrez_Gene_Id"
    } else {ID.remove <- "Hugo_Symbol"}
    
    # process
    gex.data <- gex.data %>% 
        row_to_names(row_number = 1) %>% 
        mutate_all(na_if,"") %>% 
        drop_na(!!rlang::ensym(ID.type)) %>% 
        distinct(!!rlang::ensym(ID.type),.keep_all = TRUE) %>% 
        column_to_rownames(var=ID.type) %>% 
        dplyr::select(-c(!!rlang::ensym(ID.remove))) 
    
    return(gex.data)
    
}

# define "not in" operator 
'%!in%' <- function(x,y)!('%in%'(x,y))

### functions for metagenes script (LUMA,LUMB) ###

# function to calc. the score for each metagene for each sample
mgscore <- function(metagene,metagene.def,gex.data) {
    
    # extract the gex for each metagene
    mg.ids <- filter(metagene.def, module == metagene) %>% 
        pull(entrezgene_id)
    
    # filter out the corresponding gex data
    mg.gex.data <- gex.data %>% 
        rownames_to_column(var = "entrezgene_id") %>% 
        filter(entrezgene_id %in% mg.ids) %>% 
        column_to_rownames(var = "entrezgene_id") 
    
    # calc. the score for each sample
    result <- as.data.frame(apply(mg.gex.data, 2, mean)) %>% 
        dplyr::rename(!!metagene := "apply(mg.gex.data, 2, mean)")
    
    return(result)
}

# get gex, sampleid, pam50
get_gex <- function(id,gex.data,anno) {
    # extract the gex for each gene
    single.gex <- gex.data %>% 
        select_if(~ !any(is.na(.))) %>% 
        rownames_to_column(var = "gene_id") %>% 
        filter(gene_id == id) %>% 
        column_to_rownames(var = "gene_id") %>% 
        t(.) %>% as.data.frame() %>% 
        rownames_to_column(var="sampleID")
    # combine with anno
    single.gex <- merge(single.gex, anno[,c("PAM50","sampleID")],by="sampleID") %>% relocate(PAM50,.after=sampleID)
}

# get gex, sampleid, Group for HER2p script
get_gex_hp <- function(id,gex.data,anno) {
  # extract the gex for each gene
  single.gex <- gex.data %>% 
    select_if(~ !any(is.na(.))) %>% 
    rownames_to_column(var = "gene_id") %>% 
    filter(gene_id == id) %>% 
    column_to_rownames(var = "gene_id") %>% 
    t(.) %>% as.data.frame() %>% 
    rownames_to_column(var="sampleID")
  # combine with anno
  single.gex <- merge(single.gex, anno[,c("Group","sampleID")],by="sampleID") %>% relocate(Group,.after=sampleID)
}

# summary stats.
get_stats <- function(data,grouping_var,stat_var) {
    res <- data %>% 
        group_by(!!rlang::ensym(grouping_var)) %>%
        get_summary_stats(!!rlang::ensym(stat_var), type = "mean_sd")
    return(res)
}

# ttest function that can catch errors (because of Error: data are essentially constant if I dont filter based on sd before)
e.ttest <- function(vec1,vec2,var.flag) {
    res <- tryCatch(
        expr = {
            if (var.flag) {
                t.test(vec1,vec2, var.equal = TRUE)
            } else {
                t.test(vec1,vec2, var.equal = FALSE)
            }
        },
        error = function(e) {
            list(p.value=c(NA),estimate=c(NA,NA))
        }
    )
    return(res)
}

# ttest function with equal var check
# switched out t.test() for e.ttest()
var_ttest <- function(dat1,dat2) {
    # equal var check
    vres <- var.test(dat1,dat2)
    if (!is.na(vres$p.value) & vres$p.value >= 0.05) {
        res <- e.ttest(dat1,dat2, TRUE)
    } else {
        res <- e.ttest(dat1,dat2, FALSE)
    }
    return(res)
}

# function for pairwise ttests (g1vsg2,g1vsg3)
pair_ttest <- function(data,anno=NULL,group.var,test.var,g1,g2,g3,string.output=TRUE) {
    
    # if anno data was provided it means that the data does not contain the sample pam50 annotations and is assumed to have the structure(rownames=sampleID; columns=data). The anno object does have to have the structure of sampleID as column, PAM50 (grouping var) as column
    if(!is.null(anno)) { 
        data <- merge(data %>% rownames_to_column(var="sampleID"),anno[,c("sampleID",group.var)],by="sampleID")
    }
    
    # group data
    dat1 <- data %>% 
        filter(!!rlang::ensym(group.var)==g1) %>% 
        pull(!!rlang::ensym(test.var))
    
    dat2 <- data %>% 
        filter(!!rlang::ensym(group.var)==g2) %>% 
        pull(!!rlang::ensym(test.var))
    dat3 <- data %>% 
        filter(!!rlang::ensym(group.var)==g3) %>% 
        pull(!!rlang::ensym(test.var))
    
    # test and output formatting
    if(string.output==TRUE) {
      result <- c(paste(g1,".",g2,sep=""),
                  capture.output(var_ttest(dat1,dat2)),
                  paste(g1,".",g3,sep=""),
                  capture.output(var_ttest(dat1,dat3)))
    } else {
      result <- list(
        g1.g2.res=list(var_ttest(dat1,dat2)),
        g1.g3.res=list(var_ttest(dat1,dat3)))
    }
    
    return(result)
}

# boxplot for 3 groups
three_boxplot <- function(data, group.var, test.var, g1, g2, g3, g1.col = "#d334eb", g3.col = "#34c6eb", g2.col = "#2176d5", xlab="PAM50 subtype", ylab, title, colors, ylim=NULL, subtitle=NULL) {
  
    # plot
    plot <- ggplot(data, aes(x=as.factor(.data[[group.var]]),y=.data[[test.var]],fill=as.factor(.data[[group.var]]))) +
        geom_boxplot(size=2.5, outlier.size = 7) +
        xlab(xlab) +
        ylab(ylab) +
        theme_bw() +
        theme(axis.text.x = element_text(size = 60,margin = margin(t=10)),
              axis.title.x = element_text(size = 60),
              axis.text.y = element_text(size = 55,margin = margin(r=10)),
              axis.title.y = element_text(size = 60),
              plot.title = element_text(size=50),
              legend.position = "none",
              panel.border = element_blank(), 
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black",linewidth=2),
              axis.ticks = element_line(colour = "black", linewidth = 2),
              axis.ticks.length=unit(0.5, "cm")) +
        scale_fill_manual(values=colors) +
        scale_y_continuous(breaks = scales::breaks_pretty(10)) +
        scale_x_discrete(limits = c(g2,g3,g1), labels = c(toupper(g2),toupper(g3),toupper(g1)))  # check that these are in the correct order
    if (!is.null(ylim)) {
      plot <- plot + coord_cartesian(ylim = ylim)
    } 
    
    if (!is.null(subtitle)) {
      plot <- plot + ggtitle(title,subtitle=subtitle)
    } else {
      plot <- plot + ggtitle(title)
    }
    
    return(plot)
}

# plot enrichment analysis

# labeller function for plot
set_labeller <- function(variable,value){
  return(deg.sets[value])
}

# plot enrichment analysis
pwplot <- function(res,title) {
    ggplot(res, aes(y=Gene_count,
                    x=Term,
                    fill=Adjusted.P.value)) +
    geom_bar(stat="identity") +
    coord_flip() +
    theme_bw() +
    ggtitle(title) +
    theme(text = element_text(size = 30),
          axis.text.x = element_text(margin = margin(t=5)),
          axis.text.y = element_text(margin = margin(r=5)),
          panel.border = element_blank(),          
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black",linewidth=2),
          axis.ticks = element_line(colour = "black", linewidth = 2),
          axis.ticks.length=unit(0.5, "cm"),
          legend.key.size = unit(1, 'cm')) +
          labs(fill = "Adj. p-value") +
    scale_y_continuous(limits = c(0,65),breaks = seq(0,65,5)) +
    scale_fill_gradientn(limits=c(0,0.05), colors=c("red","blue")) +
    facet_grid(rows = vars(Comp),drop=TRUE,scales="free",space="free",labeller=set_labeller)
    
}

# heatmap 

# bin continuous anno variables to color them correctly (e.g. the metagene scores)
mg_converter <- function(mg) {
    mg.class<-rep(NA,length(mg)) # global variable or what????
    mg.class[which(mg <= -2)] <- "<= -2"
    mg.class[which((mg <= -1) & (mg > -2) )] <- "-1 to -2"
    mg.class[which((mg <= -0.5) & (mg > -1) )] <- "-0.5 to -1"
    mg.class[which((mg >= -0.5) & (mg < 0.5) )] <- "-0.5 to 0.5"
    mg.class[which((mg >= 0.5) & (mg < 1) )] <- "0.5 to 1"
    mg.class[which((mg >= 1) & (mg < 2) )] <- "1 to 2"
    mg.class[which(mg >= 2)] <- ">= 2"
    return(mg.class)
}

# normality check function
select_test <- function(dat1,dat2) {
  if(diff(range(dat1)) == 0 | 
     diff(range(dat2)) == 0 ) {
    # all values identical for one or both -> cant be normally distrib.
    "non-parametric"
  } else if(shapiro.test(dat1)$p.value > 0.05 & 
            shapiro.test(dat2)$p.value > 0.05)  {
    # parametric: both normally distributed
    "parametric"
  } else {     
    # at least one group not normally distributed
    "non-parametric"
  }
}

# 2-group parametric test: t-test
var_ttest <- function(dat1,dat2) {
  # equal variance check
  levenes.res <- var.test(dat1, dat2, alternative = "two.sided")
  # t.test with or without Welchs correction
  if (!is.nan(levenes.res$p.value) & levenes.res$p.value <= 0.05) {
    t.test(dat1, dat2, var.equal = FALSE)
  } else { 
    t.test(dat1, dat2, var.equal = TRUE) 
    }
}

# 2-group non-parametric test: Mann-Whitney U
mwu_test <- function(dat1,dat2) {
  wilcox.test(dat1, dat2)
}

# combined 2-group test function
duo.test <- function(dat1,dat2) {
  test.type <- select_test(dat1,dat2)
  if (test.type=="parametric") {
    var_ttest(dat1,dat2) 
  } else if (test.type=="non-parametric") { 
    mwu_test(dat1,dat2)
  }
}
