# Function definitions for genomic analyses

################################################################################
# function statistics for mut sigs
################################################################################

# function

sig_ttest <- function(mut.data,sign){
    # vars
    hdata <- subset(mut.data, Subtype == "HER2E", select=c(sign)) 
    adata <- subset(mut.data, Subtype == "LumA", select=c(sign)) 
    bdata <- subset(mut.data, Subtype == "LumB", select=c(sign)) 
    edata <- subset(mut.data, Subtype == "ERpHER2p", select=c(sign)) 
    # for Her2 vs LumA
    # equal variance check
    if (var.test(unlist(hdata),unlist(adata), alternative = "two.sided")$p.value <= 0.05) {
        H2vsLA_ttest_result <- t.test(hdata,adata, var.equal = FALSE)
    } else {
        H2vsLA_ttest_result <- t.test(hdata,adata, var.equal = TRUE)
    }
    # save results
    print(paste("HER2E vs. LumA: ", H2vsLA_ttest_result$p.value,sep = "")) # 0.01478637
    
    # for Her2 vs LumB
    # equal variance check
    if (var.test(unlist(hdata),unlist(bdata), alternative = "two.sided")$p.value <= 0.05) {
        H2vsLB_ttest_result <- t.test(hdata,bdata, var.equal = FALSE)
    } else {
        H2vsLB_ttest_result <- t.test(hdata,bdata, var.equal = TRUE)
    }
    # save results
    print(paste("HER2E vs. LumB: ", H2vsLB_ttest_result$p.value,sep = ""))  # 0.07441933
    
    # for Her2 vs ERpHER2p
    # equal variance check
    if (var.test(unlist(hdata),unlist(edata), alternative = "two.sided")$p.value <= 0.05) {
        H2vsE_ttest_result <- t.test(hdata,edata, var.equal = FALSE)
    } else {
        H2vsE_ttest_result <- t.test(hdata,edata, var.equal = TRUE)
    }
    # save results
    print(paste("HER2E vs. ERpHER2p: ", H2vsE_ttest_result$p.value,sep = ""))  #0.01534438
    
}
