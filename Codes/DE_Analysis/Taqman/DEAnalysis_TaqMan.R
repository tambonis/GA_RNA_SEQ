################################################################################
################################################################################
# Objective: differential expression analysis using TaqMan data. 
# Author: Tiago Tambonis.
# Additional informations: 
# Date: 08/17. 
# The code is based on the scripts provided by Rapaport et al. Genome Biology 2013.
################################################################################
################################################################################

## Load TaqMan Data
load("TaqManData.Rdata")
taq.data <- as.data.frame(taq.dat, stringsAsFactors=FALSE)

## Load DE results.
load("Results_DB_2017-09-01.Rdata")

res.deseq <- results.full[["DESeq"]]
res.edger <- results.full[["edgeR"]]
res.limmaVoom <- results.full[["limmaVoom"]]
res.poseq <- results.full[["PoissonSeq"]]$res
res.bayseq <- results.full[['baySeq']]$de[,c("Likelihood", "FDR.DE")]

## List to store the data used to generate the ROC analysis.
plot.dat <- list()

###########
## DESeq
###########
## reorganize de table
rownames(res.deseq$de) <- res.deseq$de[,'id']
res.deseq$de <- res.deseq$de[,c('id', 'pval', 'padj', 'log2FoldChange', 'baseMeanA', 'baseMeanB')]

deseq.taq <- merge(res.deseq$de,
                   ##res.deseq$all.res, ## use this if using Results_tophat2.RData
                   taq.dat,
                   by.x=1, by.y='row.names')

plot.dat["DESeq"] <- list(deseq.taq)

#########
## edgeR
#########
res.edger.all <- cbind(rownames(res.edger$de$table), res.edger$de$table[,-2])
## reoder columns 
res.edger.all <- res.edger.all[,c(1,3,4,2)]
colnames(res.edger.all) <- c("ID", "Pva", "FDR", "logFC")


edger.taq <- merge(res.edger.all, taq.dat,
                   by.x='row.names', by.y='row.names')

plot.dat["edgeR"] <- list(edger.taq)

###############
## limma Voom
#############
limma.taqVoom <- merge(res.limmaVoom$tab,
                       taq.dat,
                       by.x='row.names', by.y='row.names')

plot.dat["limmaVoom"] <- list(limma.taqVoom)

############
## PoissonSeq
############
poiss.matrix <- data.frame(tt=res.poseq$tt, pval=res.poseq$pval,
                           fdr=res.poseq$fdr, logFC=res.poseq$log.fc)
rownames(poiss.matrix) <- res.poseq$gname

poseq.dat <- merge(poiss.matrix, taq.dat,
                   by.x='row.names', by.y='row.names')

plot.dat["PoissonSeq"] <- list(poseq.dat)

##############
## baySeq
##############
bayseq.taq <- merge(res.bayseq, taq.dat,
                    by.x='row.names', by.y='row.names')
## inverse M values
plot.dat["baySeq"] <- list(bayseq.taq)

##############
## Geometric approach
##############
source("GA_filter.R")
source("RPM_normalization.R")
source("Geometric_Approach.R")
load("HTSeq.RData")

group <- c(rep(1,5), rep(2,5)) #Definition of experimental conditions. 
#1 represents the condition A and 2 the condition B.

counts.dat <- HTSeq.dat
counts.dat <- matrix(as.numeric(HTSeq.dat), nrow= nrow(HTSeq.dat))
rownames(counts.dat) <- rownames(HTSeq.dat)
colnames(counts.dat) <- colnames(HTSeq.dat)
counts.dat <- counts.dat[grep('no_feature|ambiguous|too_low_aQual|not_aligned|alignment_not_unique' ,
                              rownames(counts.dat), invert=TRUE),]
counts.dat <- GA.filter(counts.dat = counts.dat)
counts.dat <- RPM_normalization(counts.dat = counts.dat)
counts.dat <- log2(counts.dat)

results <- GA(counts.dat = counts.dat, group = group) #Execute the geometric approach.
DEGA.taq <- merge(results, taq.dat, by.x='row.names', by.y='row.names')
plot.dat["Geom.Appr."] <- list(DEGA.taq) 
################################################################################
## Plots
################################################################################

#Color curves.
colr <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C")
## Number of packages analized.
kNumOfMethods <- 6

## list of adj.pval columns, FDR or relevances
qval.index <- list(DESeq=3,  edgeR=4, limmaVoom=3, PoissonSeq=4, baySeq=3, GA=2)

##############
## plot ROC
##############

kLog2Cutoffs <- c(0.5,1.0,1.5,2.0)

for (i in 1:4){
  
  kLog2Cutoff <- kLog2Cutoffs[i]
  
  plot2file=TRUE
  if(plot2file){
    setEPS()
    postscript(((paste("TaqMan_DE_analysis",kLog2Cutoff,"_", Sys.Date(),".eps", sep=''))))
  }

  PlotRocs <- function(i, dat, qval.index, logFC.index, color){
    require(pROC)
    outcome= rep(1, dim(dat)[1])
    outcome[abs(dat[,logFC.index]) <= kLog2Cutoff] =0
    if(i==1){
      roc <- plot.roc(outcome, dat[,qval.index],col=color,
                    main="ROC of TaqMan data", ylim=c(0,1.05), cex.lab=1.9,cex.main=1.85, lwd=4)
      mtext(paste("logFC cutoff= ", kLog2Cutoff, sep=''), side=3, padj=0.0, cex=2.5)
    
    }else{
      roc <- lines.roc(outcome, dat[,qval.index], add=TRUE, col=color)
    }
    return(roc)
  }


  res <- lapply(seq(kNumOfMethods), function(i) PlotRocs(i, plot.dat[[i]],
                                                       qval.index[[i]],
                                                       dim(plot.dat[[i]])[2],
                                                       colr[i]))

  names(res) <- names(plot.dat)
  legends <- lapply(seq(kNumOfMethods), function(i) paste(names(res)[i], "AUC =", format(res[[i]]$auc, digits=3), sep=' '))
  legend("bottomright", legend=legends, col=colr, lwd=5, cex=1.4, inset=c(-0,0.03))

  if(plot2file){
    dev.off()
  }
}  

#########################
## Calculate AUCs
## by changing log2 cutoff
##########################

plot2file=TRUE
if(plot2file){
  setEPS()
  postscript(paste("TaqMan_DE_analysis_AUCs", Sys.Date(), ".eps", sep=''))
}

x_AUC <- function(i, dat, qval.index, logFC.index){
  ## calculate ROC
  ## return AUC vector for a range of logFC cutoffs
  require(pROC)
  auc.res <- matrix(nrow=length(seq(0.5,2,0.1)), ncol=1)
  
  ## logFC cutoff range
  cutoff <- seq(0.5,2,0.1)
  
  for(i in seq(1:length(cutoff))){
    outcome <- rep(1, dim(dat)[1])
    outcome[abs(dat[,logFC.index]) <= cutoff[i]] =0
    
    auc.res[i] <- roc(outcome, dat[,qval.index])$auc[[1]]
  }
  return(auc.res)
}

auc.res <- sapply(seq(kNumOfMethods), function(i) x_AUC(i, plot.dat[[i]],
                                                        qval.index[[i]],
                                                        dim(plot.dat[[i]])[2])) ## TaqMan logFC is last column

colnames(auc.res) <- names(plot.dat)

## plot AUCs
plot(seq(0.5,2,0.1), auc.res[,1], type='n', main="TaqMan AUCs",
     xlab="logFC cutoff values", ylab="AUC",
     ylim=c(0.8,1), cex.lab=1.4,cex.main=1.85)

for(i in seq(dim(auc.res)[2])){
  lines(seq(0.5,2,0.1), auc.res[,i],
        lwd=4, col=colr[i])
}
legend("topleft", legend=colnames(auc.res), col=colr, lwd=5,  cex=1.1)

if(plot2file){
  dev.off()
}

#########################
## Calculate partial AUCs
## by changing log2 cutoff
##########################

plot2file=TRUE
if(plot2file){
  setEPS()
  postscript(paste("TaqMan_DE_analysis_Partial_AUCs", Sys.Date(), ".eps", sep=''))
}

x_AUC <- function(i, dat, qval.index, logFC.index){
  ## calculate ROC
  ## return AUC vector for a range of logFC cutoffs
  require(pROC)
  auc.res <- matrix(nrow=length(seq(0.5,2,0.1)), ncol=1)
  
  ## logFC cutoff range
  cutoff <- seq(0.5,2,0.1)
  
  for(i in seq(1:length(cutoff))){
    outcome <- rep(1, dim(dat)[1])
    outcome[abs(dat[,logFC.index]) <= cutoff[i]] =0
    
    auc.res[i] <- roc(outcome, dat[,qval.index], partial.auc=c(1, .8), 
                      partial.auc.focus="sp", partial.auc.correct=FALSE)$auc[[1]]
  }
  return(auc.res)
}

auc.res <- sapply(seq(kNumOfMethods), function(i) x_AUC(i, plot.dat[[i]],
                                                        qval.index[[i]],
                                                        dim(plot.dat[[i]])[2])) ## TaqMan logFC is last column

colnames(auc.res) <- names(plot.dat)

## plot AUCs
plot(seq(0.5,2,0.1), auc.res[,1], type='n', main="TaqMan partial AUCs",
     xlab="logFC cutoff values", ylab=" Partial AUC (specificity 1-0.9)",
     ylim=c(0.1,0.18), cex.lab=1.4,cex.main=1.85)

for(i in seq(dim(auc.res)[2])){
  lines(seq(0.5,2,0.1), auc.res[,i],
        lwd=4, col=colr[i])
}
legend("topleft", legend=colnames(auc.res), col=colr, lwd=5,  cex=1.1)

if(plot2file){
  dev.off()
}
