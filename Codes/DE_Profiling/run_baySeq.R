runBaySeq <- function(count.dat, conditions){
  
  if(! packageDescription("baySeq")$Version == "2.4.1"){
    stop("Wrong version of baySeq package. This script requires  baySeq-2.4.1")
  }
  
  require("baySeq")
  require("snow")
  
  ## load gene length data
  gene.length <- read.table("gene_length.txt", sep='\t', stringsAsFactors=FALSE,
                            col.names=c('gene','length', 'pos'))
  ## add spiked-in length data
  ercc.length <- read.table("ERCC92.gtf", stringsAsFactors=FALSE, sep='\t')
  ercc.length <- ercc.length[,c(1,5,9)]
  colnames(ercc.length)= colnames(gene.length)
  gene.length <- rbind(gene.length, ercc.length)
  
  ## extract seq lengths
  seqlen <- unlist(sapply(rownames(count.dat), function(x) gene.length[gene.length[,1] == x,'length'][1]))
  count.dat <- count.dat[!is.na(seqlen),]
  seqlen <- seqlen[!is.na(seqlen)]
  CD <- new("countData", data= count.dat, replicates= conditions, seglens=seqlen,
    groups = list(NDE=rep(1,length(conditions)), DE=c(rep(1,length(grep('A', conditions))),
                                                   rep(2,length(grep('B', conditions))))))

  ## estimate library scaling factor
  ## Not clear which normalization method is used
  ## assuming 0.75 quantile
  libsizes(CD) <- getLibsizes(CD)

  ## setup cluster 
  cl <- makeCluster(detectCores() - 1)
  
  CD <- getPriors.NB(CD, samplesize=1e5, estimation = 'QL', cl=cl)
  CD = getLikelihoods(CD ,pET= 'BIC', cl=cl)
  res <- topCounts(CD, group='DE', normaliseData=TRUE, likelihood=0.0)
  
  if(!is.null(cl)){
    stopCluster(cl)
  }
  return(list(dat=CD, de=res))
}

