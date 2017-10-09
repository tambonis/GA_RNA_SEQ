runLimmaQuantile <- function(count.dat, conditions,
                             q.cut=0.05, lfc=0.0, condA='condA', condB='condB'){
  
  if(! packageDescription("limma")$Version == "3.26.9"){
    stop("Wrong version of limma package. This script requires  limma-3.26.9")
  }
  
  require(limma)

  groups=as.factor(conditions)
  design <- model.matrix(~0+groups)
  colnames(design) = levels(groups)
  
  require("edgeR")
  nf <- calcNormFactors(count.dat)
  dat <- voom(count.dat, design, plot=FALSE, lib.size=colSums(count.dat) * nf)
  
  fit=lmFit(dat,design)

  contrast.matrix <- makeContrasts(condB - condA, levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  fit2 <- eBayes(fit2)
  res=decideTests(fit2,p.value=q.cut,lfc=lfc)
  tab<-topTable(fit2, adjust = "BH", number=nrow(fit2), sort.by='logFC')
  
  counts.limma=2^dat$E[rownames(tab),]
  
  val1=apply(counts.limma[,which(conditions==condA)],1,mean)
  val2=apply(counts.limma[,which(conditions==condB)],1,mean)

  tab=tab[,c("P.Value","adj.P.Val","logFC")]
  tab=as.matrix(cbind(tab,val1,val2))
  nam1=paste("mean_",condA,sep="")
  nam2=paste("mean_",condB,sep="")
  colnames(tab)[4]=nam1
  colnames(tab)[5]=nam2

  return(list(tab=tab,res=res,counts=counts.limma))

}
