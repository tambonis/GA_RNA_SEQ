runDESeq <- function(PM.dat, samples){

  if(! packageDescription("DESeq")$Version == "1.22.1"){
      stop("Wrong version of DESeq package. This script requires  DESeq-1.22.1")
    }

  require("DESeq")
  
  #Setting conditions.
  condA=c(paste("A_", seq(1,sum(as.logical(samples <= 69)))))
  condB=c(paste("B_", seq(1,sum(as.logical(samples > 69)))))
  conditions=c(rep("condA", length(condA)), rep("condB", length(condB)))
  
  kSamples <- colnames(exprs(PM.dat)[,c(samples)])
  
  targets <- data.frame(kSamples, Factor=conditions)
  
  ## Initialize new DESeq object
  cds <- newCountDataSet(exprs(PM.dat)[,c(samples)], targets$Factor)
  
  ## estimate size factor
  cds <- estimateSizeFactors(cds)
  
  ## estimate dispersion parameters (as described in paper)
  cds <- estimateDispersions(cds, method= "per-condition", fitType='local')
  
  ## differential expression
  res <- nbinomTest(cds, "condA", "condB")
  res <- res[order(res$padj), ]

  return(list(de=res))
}
