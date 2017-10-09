runedgeR <- function(PM.dat, samples){
    
    if(! packageDescription("edgeR")$Version == "3.12.1"){
      stop("Wrong version of edgeR package. This script requires edgeR-3.12.1")
    }

    require("edgeR")

    ## Initiate DGElist
    y <- DGEList(counts=exprs(PM.dat)[,c(samples)], group = PM.dat$study[c(samples)]) #DGE object.

    ## normalize
    y <- calcNormFactors(y, method="TMM")

    ## edgeR user guide suggests to set
    y <- estimateCommonDisp(y, verbose=T)
    y <- estimateTagwiseDisp(y) 
    
    ## DE test
    et <- exactTest(y, pair=c("Pickrell", "Montgomery")) 
    res <- topTags(et,n=nrow(et))

    return(list(counts=y, de=res))

}
