################################################################################
################################################################################
# Objective: calcution of performance measures using the Everaert data. 
# Author: Tiago Tambonis.
# Additional informations: 
# Date: 04/18.
# The code is based on the scripts provided by Rapaport et al. Genome Biology 2013.
################################################################################
################################################################################

library(gridExtra)

## Load DE results.
load("Results_DB_2018-02-20.Rdata")
load("HtseqCount_MAQC.RData")

dat <- read.csv("41598_2017_1617_MOESM2_ESM-5.csv", header = TRUE, row.names = 1)
names_dat <- row.names(dat)
dat <- log2(dat[,1]/dat[,2])
names(dat) <- names_dat
dat <- as.data.frame(dat)
colnames(dat) <- "log2_FC"
  
res.deseq <- results.full[["DESeq"]]
res.edger <- results.full[["edgeR"]]
res.limmaVoom <- results.full[["limmaVoom"]]
res.poseq <- results.full[["PoissonSeq"]]$res
res.bayseq <- results.full[['baySeq']]$de[,c("Likelihood", "FDR.DE")]

## List to store the data used to generate the ROC analysis.
plot.dat <- list()
all_results <- list()

###########
## DESeq
###########
## reorganize de table
rownames(res.deseq$de) <- res.deseq$de[,'id']
res.deseq$de <- res.deseq$de[,c('id', 'pval', 'padj', 'log2FoldChange', 'baseMeanA', 'baseMeanB')]

deseq.taq <- merge(res.deseq$de,
                   ##res.deseq$all.res, ## use this if using Results_tophat2.RData
                   dat,
                   by.x=1, by.y='row.names')
deseq.taq <- deseq.taq[-which(is.na(deseq.taq[,3])),]

plot.dat["DESeq"] <- list(deseq.taq)

#########
## edgeR
#########
res.edger.all <- cbind(rownames(res.edger$de$table), res.edger$de$table[,-2])
## reoder columns 
res.edger.all <- res.edger.all[,c(1,3,4,2)]
colnames(res.edger.all) <- c("ID", "Pva", "FDR", "logFC")


edger.taq <- merge(res.edger.all, dat,
                   by.x='row.names', by.y='row.names')

plot.dat["edgeR"] <- list(edger.taq)

###############
## limma Voom
#############
limma.taqVoom <- merge(res.limmaVoom$tab,
                       dat,
                       by.x='row.names', by.y='row.names')

plot.dat["limmaVoom"] <- list(limma.taqVoom)

############
## PoissonSeq
############
poiss.matrix <- data.frame(tt=res.poseq$tt, pval=res.poseq$pval,
                           fdr=res.poseq$fdr, logFC=res.poseq$log.fc)
rownames(poiss.matrix) <- res.poseq$gname

poseq.dat <- merge(poiss.matrix, dat,
                   by.x='row.names', by.y='row.names')

plot.dat["PoissonSeq"] <- list(poseq.dat)

##############
## baySeq
##############
bayseq.taq <- merge(res.bayseq, dat,
                    by.x='row.names', by.y='row.names')
## inverse M values
plot.dat["baySeq"] <- list(bayseq.taq)

##############
## Geometric approach
##############
source("GA_filter.R")
source("RPM_normalization.R")
source("Geometric_Approach.R")

group <- c(rep(1,2), rep(2,2)) #Definition of experimental conditions. 
#1 represents the condition A and 2 the condition B.

counts.dat <- maqc_data
rm(maqc_data)
  
results <- GA(counts.dat = counts.dat, group = group) #Execute the geometric approach.
results <- as.data.frame(results)
filter <- results!=0
results <- cbind(results, row.names(results))
results_order <- results[order(-results[,1])[1:sum(filter)],]
k_GA_cutoffs <- tail(results_order)[6,1]

DEGA.taq <- merge(results, dat, by.x='row.names', by.y='row.names')
plot.dat["Geom.Appr."] <- list(DEGA.taq) 

################################################################################
## Measures
################################################################################

for (z in seq(1:(length(k_GA_cutoffs)))){
  
  GA_cutoff <- k_GA_cutoffs[z]
  
  kLog2Cutoffs <- c(0.025, 0.05, 0.075, 0.1)
  
  for (j in 1:length(kLog2Cutoffs)){
    
    kLog2Cutoff <- kLog2Cutoffs[j]
    
    #Color curves.
    colr <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C")
    ## Number of packages analized.
    kNumOfMethods <- 6
    
    qval.index <- list(DESeq=3, edgeR=4, limmaVoom=3, PoissonSeq=4, baySeq=3, GA=2)
    threshold_index <- c(0.05, 0.05, 0.05, 0.05, 0.05, GA_cutoff)
    
    ## ROC
    measures <- function(i, dat, qval.index, logmix.index, color, threshold){
      
      outcome= rep(1, dim(dat)[1])
      outcome[abs(dat[,logmix.index]) <= kLog2Cutoff]=0
      if (i!=6){
        labels <- ifelse(dat[,qval.index] <= threshold, 1, 0)
        print(sum(labels))
      } else labels <- ifelse(dat[,qval.index] >= threshold, 1, 0)
      binary.labels <- labels == outcome
      
      tp <- sum((dat[,qval.index] > threshold) & binary.labels)
      fp <- sum((dat[,qval.index] > threshold) & (!binary.labels))
      tn <- sum((dat[,qval.index] <= threshold) & (!binary.labels))
      fn <-  sum((dat[,qval.index] <= threshold) & binary.labels)
      
      sensitivity <- tp/sum(binary.labels) # (True positive rate)  
      specificity <- tn/sum(!binary.labels)
      ppv <- tp/(tp+fp)
      accuracy <- (tp+tn)/(tp+fp+tn+fn)
      f1 <- 2*(ppv*sensitivity/(ppv+sensitivity))
      
      round(c(sensitivity, specificity, ppv, accuracy, f1), digits = 3)
      
    }
    
    res <- sapply(seq(kNumOfMethods), function(i) measures(i, plot.dat[[i]],
                                                           qval.index[[i]],
                                                           dim(plot.dat[[i]])[2],
                                                           colr[i], threshold_index[i]))
    res <- as.data.frame(res)
    colnames(res) <- c("DESeq", "edgeR", "limmaVoom", "PoissonSeq", "baySeq", "GA")
    rownames(res) <- c("Sensitivity.", "Specificity", "PPV", "Accuracy", "F1")
    
    all_results[as.character(kLog2Cutoff)] <- list(res)
  }
  
  for (a in seq(1:4)){
    plot2file=TRUE
    if(plot2file){
      setEPS()
      postscript(paste("Meausres_", a, "_",Sys.Date(), ".eps", sep=''))
    }
    
    res <- as.data.frame(all_results[a])
    colnames(res) <- c("DESeq", "edgeR", "limmaVoom", "PoissonSeq", "baySeq", "GA")
    grid.table(res)
    
    if(plot2file){
      dev.off()
    }
  }
}