# General functions

################################################################################
# function to process data 
################################################################################

# loads RData data file and allows to assign it directly to variable
loadRData <- function(file.path){
  load(file.path)
  get(ls()[ls() != "file.path"])
}

################################################################################
# function to cbind and fill with NA if nrow not equal
################################################################################

cbind.fill <- function(...) {
  nm <- list(...) # all arguments made into list
  nm <- lapply(nm, as.matrix) # turn into matrices, returns list object
  n <- max(sapply(nm, nrow)) # sapply return vector, n will be the highest number of rows
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

################################################################################
# "not in" operator 
################################################################################
'%!in%' <- function(x,y)!('%in%'(x,y))

################################################################################
# get stats
################################################################################

get_stats <- function(vec) {
  mean <- mean(vec, na.rm=TRUE)
  median <- median(vec, na.rm=TRUE)
  sd <- sd(vec, na.rm=TRUE)
  return(c("mean"=mean, "median"=median, "sd"=sd))
}
