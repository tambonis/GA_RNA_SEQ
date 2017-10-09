runDESeq <- function(count.dat, conditions){
  
    if(! packageDescription("DESeq")$Version == "1.22.1"){
      stop("Wrong version of DESeq package. This script requires  DESeq-1.22.1")
    }

    require("DESeq")

    kSamples <- colnames(count.dat)
    targets <- data.frame(kSamples, Factor=conditions)
    
    ## Initialize new DESeq object
    cds <- newCountDataSet(count.dat, targets$Factor)

    ## estimate size factor
    cds <- estimateSizeFactors(cds)

    ## estimate dispersion parameters (as described in paper)
    cds <- estimateDispersions(cds, method= "per-condition", fitType='local')
    
    ## differential expression
    res <- nbinomTest(cds, "condA", "condB")

    return(list(cds=cds, de=res))
}
