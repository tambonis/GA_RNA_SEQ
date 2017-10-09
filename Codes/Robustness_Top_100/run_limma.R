runLimmaQuantile <- function(PM.dat, samples,
                             q.cut=0.05, lfc=0.0){
  if(! packageDescription("limma")$Version == "3.26.9"){
    stop("Wrong version of limma package. This script requires  limma-3.26.9")
  }

  require(limma)
  require(edgeR)
  
  condA=c(paste("A_", seq(1,sum(as.logical(samples <= 69)))))
  condB=c(paste("B_", seq(1,sum(as.logical(samples > 69)))))
  conditions=c(rep("condA", length(condA)), rep("condB", length(condB)))
  
  groups=as.factor(conditions)
  design <- model.matrix(~0+groups)
  colnames(design) = levels(groups)
  
  nf <- calcNormFactors(exprs(PM.dat)[,c(samples)])
  dat <- voom(exprs(PM.dat)[,c(samples)], design, plot=FALSE, 
              lib.size=colSums(exprs(PM.dat)[,c(samples)]) * nf)
  
  fit=lmFit(dat,design)
  
  contrast.matrix <- makeContrasts("condB - condA", levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  fit2 <- eBayes(fit2)
  res=decideTests(fit2,p.value=q.cut,lfc=lfc)
  tab<-topTable(fit2, adjust = "BH", number=nrow(fit2), sort.by='logFC')
  
  tab <- tab[order(tab$adj.P.Val), ]

  return(list(tab=tab))

}
