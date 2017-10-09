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
load("Results_DB_2017-09-01.Rdata")
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
## PoissonSeq
##############
res.poseq <- results.full[["PoissonSeq"]]$res
poiss.matrix <- data.frame(tt=res.poseq$tt, pval=res.poseq$pval,
                           fdr=res.poseq$fdr, logFC=res.poseq$log.fc)
row.names(poiss.matrix) <- res.poseq$gname

PoissonSeq <- as.data.frame(poiss.matrix$pval)
row.names(PoissonSeq) <- rownames(poiss.matrix)

PoissonSeq <- unlist(-log10(PoissonSeq))
PoissonSeq[is.infinite(PoissonSeq)] <- max(PoissonSeq[is.finite(PoissonSeq)]) +1
PoissonSeq <- as.data.frame(PoissonSeq)
row.names(PoissonSeq) <- rownames(poiss.matrix)

Poisson_table <- merge(FC, PoissonSeq, by='row.names')
Poisson_table <- Poisson_table[with(Poisson_table, order(-PoissonSeq)), ]
colnames(Poisson_table) <- c("Names", "log2FC", "PoissonSeq_Pvalues")

row.names(Poisson_table) <- Poisson_table$Names
Poisson_table<- Poisson_table[,-1]

##################
#Color definition 
##################
cols <- rep("gray", dim(Poisson_table)[1])
cols[1:100] <- "red"

sizes <- rep(1, dim(Poisson_table)[1])
sizes[1:100] <- 4

Poisson_table <- cbind(Poisson_table, cols)
Poisson_table$cols <- as.character(Poisson_table$cols)
Poisson_table <- cbind(Poisson_table, sizes)

plot2file=TRUE
if(plot2file){
  #setEPS()
  #postscript(paste("PoissonSeq_fold_change_", Sys.Date(), ".eps", sep=''), width = 10, height = 10)
  jpeg(paste("PoissonSeq_fold_change_", Sys.Date(), ".jpg", sep=''))
}

plot(Poisson_table[,1], Poisson_table[,2], main="",
     xlab="Fold-change", ylab="-log10(p)", cex = as.numeric(Poisson_table$sizes), 
     pch=20, col=Poisson_table$cols, cex.lab=1.5,cex.main=1.85, lwd=4) 

legends <- c("Top 100")
legend("bottomleft", legend=legends, col=c("red"), lwd=5,  cex=1.2, xpd = TRUE, inset = c(0.1,0))

if(plot2file){
  dev.off()
}
