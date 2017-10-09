runPoissonSeq <- function(count.dat, conditions, count.cut){
  
  if(! packageDescription("PoissonSeq")$Version == "1.1.2"){
    stop("Wrong version of PoissonSeq package. This script requires  PoissonSeq-1.1.2")
  }
  
  require("PoissonSeq")

  y=conditions
  y[y=='condA'] =1
  y[y=='condB'] =2
  y <- as.numeric(y)

  condition.type <- 'twoclass'
  
  dat <- list(n=count.dat,
              y=y,
              type=condition.type,
              pair=FALSE,
              gname=rownames(count.dat))
  para <- list(ct.sum=count.cut, ct.mean=count.cut/2,npermu=500)
  
  res <- PS.Main(dat=dat, para=para)
  lib.factors <- PS.Est.Depth(count.dat)
  return(list(res=res,norm.factors=lib.factors ))
}
