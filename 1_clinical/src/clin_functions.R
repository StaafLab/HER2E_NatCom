# Function definitions for clinical data analyses

################################################################################
# functions
################################################################################
# define "not in" operator 
'%!in%' <- function(x,y)!('%in%'(x,y))

# loads RData data file and allows to assign it directly to variable
loadRData <- function(file.path){
  load(file.path)
  get(ls()[ls() != "file.path"])
}

# function to ttest two groups
TwoGroups.ttest <- function(g1.dat,g2.dat) {
  if (var.test(unlist(g1.dat),unlist(g2.dat), alternative = "two.sided")$p.value <= 0.05) {
    
    result <- t.test(g1.dat, g2.dat, var.equal = FALSE)
  } else {
    result <- t.test(g1.dat, g2.dat, var.equal = TRUE)  
  }
  
  return(result)
}

# function for 2x2 ct testing
test.2x2ct <- function(data,var.df,var) {
  
  # for Her2 vs LumA
  luma.ct <- table(data$PAM50,data[[var]])[c("Her2","LumA"),]
  var.df$LUMA.pval <- fisher.test(luma.ct)$p.value
  
  # for Her2 vs LumB
  lumb.ct <- table(data$PAM50,data[[var]])[c("Her2","LumB"),]
  var.df$LUMB.pval <- fisher.test(lumb.ct)$p.value
  
  groups <- c("Her2","LumA","LumB")
  cols <- c("HER2E","LUMA","LUMB")
  
  for (i in 1:3) {
    type <- groups[i]
    col <- cols[i]
    type.dat <- data[which(data$PAM50==type),]
    type.dat$PAM50 <- droplevels(type.dat$PAM50)
    type.counts <- table(type.dat$PAM50,type.dat[[var]])
    
    # count column
    if (col=="HER2E") {
      var.df[[paste(col,"(ref)",sep="")]][1] <- type.counts[1]
      var.df[[paste(col,"(ref)",sep="")]][2] <- type.counts[2]
    } else {
      var.df[[col]][1] <- type.counts[1]
      var.df[[col]][2] <- type.counts[2]
    }
    # % column
    var.df[[paste(col,".%",sep="")]][1] <- round(type.counts[1]/sum(type.counts)*100)
    var.df[[paste(col,".%",sep="")]][2] <- round(type.counts[2]/sum(type.counts)*100)
  } 
  return(var.df)
} 

# function for cont var testing (size and age)
test.cont <- function(data,var.df,var) {
  
  # for Her2 vs LumA
  var.df$LUMA.pval <-     TwoGroups.ttest(data[Her2.ids,var],data[LumA.ids,var])$p.value
  
  # for Her2 vs LumB
  var.df$LUMB.pval <-     TwoGroups.ttest(data[Her2.ids,var],data[LumB.ids,var])$p.value
  
  # descriptive stats
  # mean
  var.df$`HER2E(ref)` <-  mean(data[which(data$PAM50=="Her2" & !is.na(data[[var]])),][[var]])
  var.df$LUMA <- mean(data[which(data$PAM50=="LumA" & !is.na(data[[var]])),][[var]])
  var.df$LUMB <- mean(data[which(data$PAM50=="LumB" & !is.na(data[[var]])),][[var]])
  # sd
  var.df$`HER2E.%` <-  sd(data[which(data$PAM50=="Her2" & !is.na(data[[var]])),][[var]])
  var.df$`LUMA.%` <- sd(data[which(data$PAM50=="LumA" & !is.na(data[[var]])),][[var]])
  var.df$`LUMB.%` <- sd(data[which(data$PAM50=="LumB" & !is.na(data[[var]])),][[var]])
  
  return(var.df)
}

# function for 2x3 ct testing
test.2x3ct <- function(data,var.df,var) { 

  # for Her2 vs LumA
  luma.ct <- table(data$PAM50,data[[var]])[c("Her2","LumA"),]
  var.df$LUMA.pval <- fisher.test(luma.ct)$p.value
  
  # for Her2 vs LumB
  lumb.ct <- table(data$PAM50,data[[var]])[c("Her2","LumB"),]
  var.df$LUMB.pval <- fisher.test(lumb.ct)$p.value
  
  groups <- c("Her2","LumA","LumB")
  cols <- c("HER2E","LUMA","LUMB")
  
  for (i in 1:3) {
    type <- groups[i]
    col <- cols[i]
    type.dat <- data[which(data$PAM50==type),]
    type.dat$PAM50 <- droplevels(type.dat$PAM50)
    type.counts <- table(type.dat$PAM50,type.dat[[var]])
    
    # count column
    if (col=="HER2E") {
      var.df[[paste(col,"(ref)",sep="")]][1] <- type.counts[1]
      var.df[[paste(col,"(ref)",sep="")]][2] <- type.counts[2]
      var.df[[paste(col,"(ref)",sep="")]][3] <- type.counts[3]
    } else {
      var.df[[col]][1] <- type.counts[1]
      var.df[[col]][2] <- type.counts[2]
      var.df[[col]][3] <- type.counts[3]
    }
    # % column
    var.df[[paste(col,".%",sep="")]][1] <- round(type.counts[1]/sum(type.counts)*100)
    var.df[[paste(col,".%",sep="")]][2] <- round(type.counts[2]/sum(type.counts)*100)
    var.df[[paste(col,".%",sep="")]][3] <- round(type.counts[3]/sum(type.counts)*100)
  }
  
  return(var.df) 
}

# function that creates KM plot for specified OM for 3 Her2e, luma , lumb
KMplot <- function(OM,OMbin,OMstring,group.cohort.version,sdata) {
    
    # surv object
    data.surv <- Surv(OM, OMbin) 
    
    # fit
    fit <- survminer::surv_fit(data.surv~PAM50, data=sdata, conf.type="log-log") 
    print(head(fit))
    # weird bug: survival::survfit() cant be passed data in function call ?! so i use survminer::surv_fit()
    #survdiff(data.surv ~ PAM50, data = sdata) 
    
    plot <- ggsurvplot(
        fit,
        censor.size = 8,
        censor.shape = "|",
        size = 5,
        risk.table = FALSE,       
        pval = TRUE,
        pval.size = 8,
        pval.coord = c(0,0.1),
        conf.int = FALSE,         
        xlim = c(0,max(OM[is.finite(OM)])),         
        xlab = paste(OMstring," (years)", sep = ""),
        ylab = paste(OMstring," event probability", sep = ""), # ggf just label as "event probability"
        ylim = c(0,1),
        palette = c("#2176d5", "#34c6eb", "#d334eb"), 
        legend = c(0.9,0.96),
        ggtheme = theme(legend.title = element_text(size=25), #20
                        legend.key.size = unit(0.5,"cm"), 
                        legend.text = element_text(size = 25), #20
                        axis.text.x = element_text(size = 25), #20
                        axis.title.x = element_text(size = 30), #25
                        axis.text.y = element_text(size = 25), #20
                        axis.title.y = element_text(size = 30),
                        plot.title = element_text(size=30),
                        panel.background = element_blank(),
                        panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),
                        axis.line = element_line(colour = "black",linewidth=2),
                        axis.ticks = element_line(colour = "black", linewidth = 2),
                        axis.ticks.length=unit(0.25, "cm")),
        title= paste(OMstring, ": ",group.cohort.version, sep = ""),
        legend.title = "Subtypes",
        legend.labs = c(paste("LUMA"," (",table(sdata[!is.na(OM),]$PAM50)[1],")",sep = ""),
                        paste("LUMB"," (",table(sdata[!is.na(OM),]$PAM50)[2],")",sep = ""),
                        paste("HER2E"," (",table(sdata[!is.na(OM),]$PAM50)[3],")",sep = "")),
        break.x.by = 1, # break X axis in time intervals of x (what is nicest here? maybe 365)
        break.y.by = 0.1)
    
    return(plot)
    
}

# function that creates KM plot for specified OM for two groups
TwoGroup.KMplot <- function(OM,OMbin,OMstring,group.cohort.version,sdata,comp.var,palette,legend.labs) {
  
  # surv object
  data.surv <- Surv(OM, OMbin) 
  
  # fit
  fit <- survminer::surv_fit(data.surv~get(comp.var), data=sdata, conf.type="log-log") # weird bug: survival::survfit() cant be passed data in function call ?! so i use survminer::surv_fit()
  #survdiff(data.surv ~ PAM50, data = sdata) 
  
  plot <- ggsurvplot(
    fit,
    censor.size = 8,
    censor.shape = "|",
    size = 5,
    risk.table = FALSE,       
    pval = TRUE,
    pval.size = 8,
    pval.coord = c(0,0.1),
    conf.int = FALSE,         
    xlim = c(0,max(OM[is.finite(OM)])),         
    xlab = paste(OMstring," (years)", sep = ""),
    ylab = paste(OMstring," event probability", sep = ""), # ggf just label as "event probability"
    ylim = c(0,1),
    palette = palette, 
    legend = c(0.9,0.96),
    ggtheme = theme(legend.title = element_text(size=25), #20
                    legend.key.size = unit(0.5,"cm"), 
                    legend.text = element_text(size = 25), #20
                    axis.text.x = element_text(size = 25), #20
                    axis.title.x = element_text(size = 30), #25
                    axis.text.y = element_text(size = 25), #20
                    axis.title.y = element_text(size = 30),
                    plot.title = element_text(size=30),
                    panel.background = element_blank(),
                    panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(colour = "black",linewidth=2),
                    axis.ticks = element_line(colour = "black", linewidth = 2),
                    axis.ticks.length=unit(0.25, "cm")),
    title= paste(OMstring, ": ",group.cohort.version, sep = ""),
    legend.title = "Subtypes",
    legend.labs = legend.labs,
    break.x.by = 1, # break X axis in time intervals of x (what is nicest here? maybe 365)
    break.y.by = 0.1)
  
  return(plot)
  
}

# function that created univariate Cox proportional hazards model forest plot
unicox <- function(data,surv,title) {
    
    # Model construction
    main.pam50 <- coxph(surv~PAM50, data=data)
    
    # Result
    res <- summary(main.pam50)
    #result <- ggplotify::as.ggplot(arrangeGrob(textGrob(res$call),tableGrob(res$coefficients),tableGrob(res$conf.int)))

    # forest 
    plot <- ggforest(main.pam50,fontsize = 2,data=data,main=title)
    
    return(out <- list("plot" = plot, "result" = res))
}

# function that created multivariate Cox proportional hazards model forest plot
mvcox <- function(data,surv,title) {
    
    # Model construction 
    # parameters to incl: PAM50, Age, Grade, TumSize
    main.all <- coxph(surv~PAM50+Age+LN+TumSize+Grade, data=data) 
    
    # Result
    res <- summary(main.all)
    #result <- ggplotify::as.ggplot(arrangeGrob(textGrob(res$call),tableGrob(res$coefficients),tableGrob(res$conf.int)))
    
    # Plot forest 
    plot <- ggforest(main.all,
                     fontsize = 2,
                     cpositions = c(0.01,0.13,0.35),
                     data=data,
                     main=title) 
    
    return(out <- list("plot" = plot, "result" = res))
}

# function that creates KM plot for HER2p for specified OM
HER2p_KMplot <- function(OM,OMbin,OMstring,group.cohort.version,sdata) {
    
    # surv object
    data.surv <- Surv(OM, OMbin) 
    
    # fit
    fit <- survminer::surv_fit(data.surv~PAM50, data=sdata, conf.type="log-log") # weird bug: survival::survfit() cant be passed data in function call ?! so i use survminer::surv_fit()
    #survdiff(data.surv ~ PAM50, data = sdata) 
    
    plot <- ggsurvplot(
        fit,
        censor.size = 8,
        censor.shape = "|",
        size = 5,
        risk.table = FALSE,       
        pval = TRUE,
        pval.size = 8,
        pval.coord = c(0,0.1),
        conf.int = FALSE,         
        xlim = c(0,max(OM[is.finite(OM)])),         
        xlab = paste(OMstring," (days)", sep = ""),
        ylab = paste(OMstring," event probability", sep = ""), # ggf just label as "event probability"
        ylim = c(0,1),
        palette = c("#d334eb", "#2176d5"), 
        legend = c(0.9,0.96),
        ggtheme = theme(legend.title = element_text(size=25), #20
                        legend.key.size = unit(0.5,"cm"), 
                        legend.text = element_text(size = 25), #20
                        axis.text.x = element_text(size = 25), #20
                        axis.title.x = element_text(size = 30), #25
                        axis.text.y = element_text(size = 25), #20
                        axis.title.y = element_text(size = 30),
                        plot.title = element_text(size=30),
                        panel.background = element_blank(),
                        panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),
                        axis.line = element_line(colour = "black",linewidth=2),
                        axis.ticks = element_line(colour = "black", linewidth = 2),
                        axis.ticks.length=unit(0.25, "cm")),
        title= paste(OMstring, ": ",group.cohort.version, sep = ""),
        legend.title = "Subtypes",
        legend.labs = c(paste("HER2E"," (",table(sdata[!is.na(OM),]$PAM50)[1],")",sep = ""),
                        paste("nonHER2E"," (",table(sdata[!is.na(OM),]$PAM50)[2],")",sep = "")),
        break.x.by = 500, # break X axis in time intervals of x (what is nicest here? maybe 365)
        break.y.by = 0.1)
    
    print(plot)
    }