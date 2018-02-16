################################################################################
################################################################################
# Objective: differential expression analysis using the ERCC data.
# Author: Tiago Tambonis and Marcelo Boareto.
# Additional informations: 
# Date: 08/17. 
# The code is based on the scripts provided by Rapaport et al. Genome Biology 2013.
################################################################################
################################################################################

## Load ERCC data
ercc <- read.table("ERCC_Controls_Analysis.txt", stringsAsFactors=FALSE, 
                   header=T, sep="\t", row.names=2)

## Load DE results.
load("Results_DB_2017-09-01.Rdata")
res.deseq <- results.full[["DESeq"]]
res.edger <- results.full[["edgeR"]]
res.limmaVoom <- results.full[["limmaVoom"]]
res.poseq <- results.full[["PoissonSeq"]]$res
res.bayseq <- results.full[['baySeq']]$de[,c("Likelihood", "FDR.DE")]

## List to store the data used to generate the ROC analysis.
plot.dat <- list()

############
## DESeq
###########
rownames(res.deseq$de) <- res.deseq$de[,'id']
res.deseq$de <- res.deseq$de[,c('id', 'pval', 'padj', 'log2FoldChange', 'baseMeanA', 'baseMeanB')]

deseq.taq <- merge(res.deseq$de,
                   ## res.deseq$all.res,
                   ercc,
                   by.x='row.names', by.y='row.names')

plot.dat["DESeq"] <- list(deseq.taq)

#########
## edgeR
#########

res.edger.all <- cbind(rownames(res.edger$de$table), res.edger$de$table[,-2])
## reoder columns 
res.edger.all <- res.edger.all[,c(1,3,4,2)]
colnames(res.edger.all) <- c("ID", "Pva", "FDR", "logFC")

edger.taq <- merge(res.edger.all, ercc,
                   by.x='row.names', by.y='row.names')

plot.dat["edgeR"] <- list(edger.taq)

###############
## limma
#############
limma.taqVoom <- merge(res.limmaVoom$tab,
                       ercc,
                       by.x='row.names', by.y='row.names')

plot.dat["limmaVoom"] <- list(limma.taqVoom)

############
## PoissonSeq
############
poiss.matrix <- data.frame(tt=res.poseq$tt, pval=res.poseq$pval,
                           fdr=res.poseq$fdr, logFC=res.poseq$log.fc)
rownames(poiss.matrix) <- res.poseq$gname

poseq.dat <- merge(poiss.matrix, ercc,
                   by.x='row.names', by.y='row.names')

plot.dat["PoissonSeq"] <- list(poseq.dat)

##############
## baySeq
##############
bayseq.dat <- merge(res.bayseq, ercc,
                    by.x='row.names', by.y='row.names')

plot.dat['baySeq'] <- list(bayseq.dat)

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
DEGA.taq <- merge(results, ercc, by.x='row.names', by.y='row.names')
plot.dat["Geom.Appr."] <- list(DEGA.taq) 

################################################################################
## Plots
################################################################################

plot2file=TRUE
if(plot2file){
  setEPS()
  postscript(paste("ERCC_DEAnalysis", Sys.Date(), ".eps", sep=''))
}

#Color curves.
colr <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C")
## Number of packages analized.
kNumOfMethods <- 6

qval.index <- list(DESeq=4, edgeR=4, limmaVoom=3, PoissonSeq=4, baySeq=2, GA=2)

## ROC
require(pROC)
PlotRocs <- function(i, dat, qval.index, logmix.index, color){
  outcome= rep(1, dim(dat)[1])
  outcome[dat[,logmix.index] == 0] =0
  if(i==1){
    roc <- plot.roc(outcome, dat[,qval.index],col=color,main="ROC of ERCC spike-in data", 
                    cex.lab=2.3,cex.main=1.85, lwd=4, cex=2.5)
    
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
legends <- lapply(seq(kNumOfMethods), 
                  function(i) paste(names(res)[i], "AUC =", format(res[[i]]$auc, 
                                                                   digits=3), sep=' '))
legend("bottomright", legend=legends, col=colr, lwd=5,  cex=1.2)

if(plot2file){
  dev.off()
}

##############
# Partial AUC
##############

#Color curves.
colr <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C")
## Number of packages analized.
kNumOfMethods <- 6

qval.index <- list(DESeq=4, edgeR=4, limmaVoom=3, PoissonSeq=4, baySeq=2, GA=2)

## ROC
require(pROC)
PlotRocs <- function(i, dat, qval.index, logmix.index, color){
  outcome= rep(1, dim(dat)[1])
  outcome[dat[,logmix.index] == 0] =0
  if(i==1){
    roc <- roc(outcome, dat[,qval.index], partial.auc=c(1, .9), 
               partial.auc.focus="sp", partial.auc.correct=FALSE)$auc[[1]]

  }else{
    roc <- roc(outcome, dat[,qval.index],  partial.auc=c(1, .9), 
                     partial.auc.focus="sp", partial.auc.correct=FALSE)$auc[[1]]
  }
  return(roc)
}

res <- lapply(seq(kNumOfMethods), function(i) PlotRocs(i, plot.dat[[i]],
                                                       qval.index[[i]],
                                                       dim(plot.dat[[i]])[2],
                                                       colr[i]))

plot2file=TRUE
if(plot2file){
  setEPS()
  postscript(paste("AUC_partial_ERCC", Sys.Date(), ".eps", sep=''))
}

library(gridExtra)
res <- round(as.numeric(res), digits = 3)
names(res) <- names(plot.dat)
res <- as.data.frame(res)
colnames(res) <- "AUC partial"
grid.table(res)

if(plot2file){
  dev.off()
}
