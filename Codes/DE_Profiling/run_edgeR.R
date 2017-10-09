runedgeR <- function(count.dat, conditions){
  
    if(! packageDescription("edgeR")$Version == "3.12.1"){
      stop("Wrong version of edgeR package. This script requires edgeR-3.12.1")
    }
    
    require("edgeR")
  
    ## Initiate DGElist
    y <- DGEList(counts=count.dat, group=factor(conditions))

    ## normalize
    y <- calcNormFactors(y)

    ## edgeR user guide suggests to set
    ## dispersion to 0.4 when there are no replicates
    if(2 == length(conditions)){
      y$common.dispersion = 0.4

    }else{
      ## estimate common dispersion
      y <- estimateCommonDisp(y)

      ## estimate gene specific dispersion
      y <- estimateTagwiseDisp(y)
    }
    ## DE test
    et <- exactTest(y)
    res <- topTags(et,n=nrow(et))

    return(list(counts=y, de=res))

}
