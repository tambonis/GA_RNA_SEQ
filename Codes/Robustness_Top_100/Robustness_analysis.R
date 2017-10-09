################################################################################
################################################################################
# Objective: Robustness analysis.
# Author: Tiago Tambonis.
# Additional informations: 
# Date: 08/17. 
################################################################################
################################################################################

library("Hmisc")

load("DESeq_overlaps.RData")
load("edgeR_overlaps.RData")
load("GA_overlaps.RData")
load("PoissonSeq_overlaps.RData")
load('limma_overlaps.RData')

colr <- c("#990000", "#000000", "#0000ff", "#006600", "#ffff00")

plot2file=TRUE
if(plot2file){
  setEPS()
  postscript(paste("Robustness_Analysis", Sys.Date(), ".eps", sep=''))
}

par(fg = "black")
errbar(1:8, rep(1000,8), rep(1000,8) + 1, rep(1000,8) - 1, 
       main="Robustness study", type = 'b', col="black", ylim=c(0,105), xlab='Sample Size',
       ylab='Average overlap', cex.lab=1.9, xaxt = 'n')
axis(1, at=1:8, labels=c(15,20,25,30,35,40,45,50))
title("Robustness analysis", sub = "", cex.main = 2, font.main= 2)
legends <- c("Geom. appr.", "edgeR", "DESeq", "PoissonSeq", "limma")
legend("bottomright", legend=legends, col=colr, lwd=4.5, cex=1.0, inset=c(-0,0.03))
par(fg = colr[1])
means <- apply(k.overlap.GA, 1, mean)
std <- apply(k.overlap.GA, 1, function(x) sqrt(var(x)))
errbar(1:8, means, means + std, means - std, type = 'b', col=colr[1], ylim=c(0,110), add=TRUE, lwd=2)

par(fg = colr[2])
means <- apply( k.overlap.edger, 1, mean)
std <- apply( k.overlap.edger, 1, function(x) sqrt(var(x)))
errbar(1:8, means, means + std, means - std, type = 'b', col=colr[2], ylim=c(0,100), add=TRUE, lwd=2)

par(fg = colr[3])
means <- apply(k.overlap.deseq, 1, mean)
std <- apply(k.overlap.deseq, 1, function(x) sqrt(var(x)))
errbar(1:8, means, means + std, means - std, type = 'b', col=colr[3], ylim=c(0,100), add=TRUE, lwd=2)

par(fg = colr[4])
means <- apply(k.overlap.poissonseq, 1, mean)
std <- apply(k.overlap.poissonseq, 1, function(x) sqrt(var(x)))
errbar(1:8, means, means + std, means - std, type = 'b', col=colr[4], ylim=c(0,100), add=TRUE, lwd=2)

par(fg = colr[5])
means <- apply( k.overlap.limma, 1, mean)
std <- apply( k.overlap.limma, 1, function(x) sqrt(var(x)))
errbar(1:8, means, means + std, means - std, type = 'b', col=colr[5], ylim=c(0,100), add=TRUE, lwd=2)

if(plot2file){
  dev.off()
}
