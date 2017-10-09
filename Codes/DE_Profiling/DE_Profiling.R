################################################################################
################################################################################
#Objective: DE analysis.
#Author: Tiago Tambonis.
#Additional informations: 
#Date: 08/17. 
################################################################################
################################################################################

#################
##Load 
#################
source("run_edgeR.R")
source("run_baySeq.R")
source("run_PoisSeq.R")
source("run_DESeq.R")
source("run_limma.R")

load("HTSeq.RData")
counts.dat <- HTSeq.dat
counts.dat <- matrix(as.numeric(HTSeq.dat), nrow= nrow(HTSeq.dat))
rownames(counts.dat) <- rownames(HTSeq.dat)
colnames(counts.dat) <- colnames(HTSeq.dat)
## remove non-gene entries
## see http://www-huber.embl.de/users/anders/HTSeq/doc/count.html
## for the 'special counters' that are added at the end
counts.dat <- counts.dat[grep('no_feature|ambiguous|too_low_aQual|not_aligned|alignment_not_unique' ,
                              rownames(counts.dat), invert=TRUE),]

results.full <- list() #List that will contain all results.

#################
##edgeR 
#################
results.full["edgeR"] <- list(runedgeR(counts.dat, c(rep("condA", 5),
                                                     rep("condB", 5))))

##############
##baySeq 
##############
results.full["baySeq"] <- list(runBaySeq(counts.dat, c(rep("condA", 5),
                                                       rep("condB", 5))))

##############
##PoissonSeq
##############
count.cut=0
results.full["PoissonSeq"] <- 
  list(runPoissonSeq(counts.dat,c(rep("condA", 5), rep("condB", 5)), count.cut))

#################
##DESeq 
#################
results.full["DESeq"] <- list(runDESeq(counts.dat, c(rep("condA", 5),
                                                     rep("condB", 5))))

##################
##limma+voom
##################
results.full["limmaVoom"] <- list(runLimmaQuantile(counts.dat, c(rep("condA", 5), rep("condB", 5))))

##############
##Salve results
##############
Script.Name <- "DE_Profilling.R"

save(results.full, file=paste("Results_DB", "_",Sys.Date(),".Rdata", sep=''))

print(paste(Script.Name, ", DONE.", sep=""))