################################################################################
################################################################################
# Objective: Fold-change study.
# Author: Tiago Tambonis.
# Additional informations: 
# Date: 08/17. 
################################################################################
################################################################################

##################
#Load
##################
source("GA_filter.R")
source("RPM_normalization.R")
source("Geometric_Approach.R")
load("HTSeq.RData")

counts.dat <- HTSeq.dat
counts.dat <- matrix(as.numeric(HTSeq.dat), nrow= nrow(HTSeq.dat))
rownames(counts.dat) <- rownames(HTSeq.dat)
colnames(counts.dat) <- colnames(HTSeq.dat)
counts.dat <- counts.dat[grep('no_feature|ambiguous|too_low_aQual|not_aligned|alignment_not_unique' ,
                              rownames(counts.dat), invert=TRUE),]

##################
#Fold-change calculation
##################
mean_a <- apply(counts.dat[,1:5], 1, mean)
mean_b <- apply(counts.dat[,6:10], 1, mean)

FC <- log2(mean_a) - log2(mean_b)

##############
## Geometric Approach
##############
counts.dat <- GA.filter(counts.dat = counts.dat)
counts.dat <- RPM_normalization(counts.dat = counts.dat)
counts.dat <- log2(counts.dat)
group <- c(rep(1,5), rep(2,5))
results.GA <- GA(counts.dat = counts.dat, group = group)

GA_table <- merge(FC, results.GA, by='row.names')
GA_table <- GA_table[with(GA_table, order(-Results)), ]

row.names(GA_table) <- GA_table$Row.names
GA_table <- GA_table[,-1]

##################
#Color definition 
##################
cols <- rep("gray", dim(GA_table)[1])
cols[1:100] <- "red"

sizes <- rep(1, dim(GA_table)[1])
sizes[1:100] <- 4

GA_table <- cbind(GA_table, cols)
GA_table$cols <- as.character(GA_table$cols)
GA_table <- cbind(GA_table, sizes)

plot2file=TRUE
if(plot2file){
  #setEPS()
  #postscript(paste("GA_fold_change_", Sys.Date(), ".eps", sep=''), width = 10, height = 10)
  jpeg(paste("GA_fold_change_", Sys.Date(), ".jpg", sep=''))
}

plot(GA_table[,1], GA_table[,2], main="",
     xlab="Fold-change", ylab="Relevance", cex = as.numeric(GA_table$sizes), 
     pch=20, col=GA_table$cols, cex.lab=1.5,cex.main=1.85, lwd=4) 
lines(GA_table[,1], rep(1,dim(GA_table)[1]), type="l") 

legends <- c("Top 100")
legend("top", legend=legends, col=c("red"), lwd=5,  cex=1.2)

if(plot2file){
  dev.off()
}
