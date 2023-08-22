rm.outlier <- function(x,percent){
  tmp <- (x>=quantile(x,percent))
  Q13 <- quantile(x[tmp],c(0.25,0.75))
  if (Q13[2]==Q13[1]) Q13[2] <- Q13[1]+1
  upper <- max(Q13[2]*3-2*Q13[1],1)
  lower <- max(Q13[1]*3-2*Q13[2],0)
  outlier.index <- tmp & (x>upper | x<lower)
  mean.noout <- mean(x[tmp & (!outlier.index)])
  x[outlier.index] <- round(mean.noout)
  cat(sum(outlier.index),"\n")
  x
}

sub.sample <- function(x){
  library.size <- x[1]
  x <- x[-1]
  sample.count <- function(a,library.size){
    rpois(1,a/library.size)
  }
  apply(matrix(x), 1, sample.count,library.size = library.size)
}

my.min <- function(x){
  if (sum(is.na(x))==length(x)){
    NA
  } else {
    min(x,na.rm = T)
  }
}

my.max <- function(x){
  if (sum(is.na(x))==length(x)){
    NA
  } else {
    max(x,na.rm = T)
  }
}

`%or%` = function(a, b) {
  cmp = function(a,b) if (identical(a, FALSE) ||
                          is.null(a) ||
                          is.na(a) ||
                          is.nan(a) ||
                          length(a) == 0) b else a
  
  if (length(a) > 1)
    mapply(cmp, a, b)
  else
    cmp(a, b)
}