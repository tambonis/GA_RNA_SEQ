runPoissonSeq <- function(PM.dat, samples){

  if(! packageDescription("PoissonSeq")$Version == "1.1.2"){
    stop("Wrong version of PoissonSeq package. This script requires  PoissonSeq-1.1.2")
  }

  require("PoissonSeq")

  count.cut=0
  
  condA=c(paste("A_", seq(1,sum(as.logical(samples <= 69)))))
  condB=c(paste("B_", seq(1,sum(as.logical(samples > 69)))))
  conditions=c(rep("condA", length(condA)), rep("condB", length(condB)))
  
  y=conditions
  y[y=='condA'] =1
  y[y=='condB'] =2
  y <- as.numeric(y)
  
  condition.type <- 'twoclass'
  
  dat <- list(n=exprs(PM.dat)[,c(samples)],
              y=y,
              type=condition.type,
              pair=FALSE,
              gname=rownames(exprs(PM.dat)[,c(samples)]))
  para <- list(ct.sum=0, pow.file="")
  
  res <- PS.Main(dat=dat, para=para)
  
  res$rnk <- sign(res$log.fc) * res$tt
  rnk <- res[, c("rnk"), drop=FALSE]
  
  return(list(res=rnk))
  
}
