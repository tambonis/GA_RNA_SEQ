################################################################################
################################################################################
# Última modificação: 27/10/16
# Objetivo: analyze the feature selection associated to Suvrel Method.
# Autor: Tiago Tambonis.
################################################################################
################################################################################

############
## Load.
###########
source("Suvrel_normalization.R") # Load the Mean 0, variance 1 normalization function.
source("Suvrel_cals.R") #Load the Suvrel implementation function.
library(tweeDEseqCountData) 
library(edgeR)

#RNA-seq data
data(pickrell) #pickrell.eset of Lymphoblastoid RNA-Seq, African (Nigeria) pop. of 69.
data(montgomery) #montgomery.eset of Lymphoblastoid RNA-Seq, Euro (Utah) pop. of 60.
PM.eset <- combine(pickrell.eset, montgomery.eset) #combined dataset, to test normalization methods.
countsPM <- exprs(PM.eset) #will test for DE between the 2 populations.
d.PM <- DGEList(counts=countsPM, group = PM.eset$study)
keep <- rowSums(cpm(d.PM)>0) >= 1 #Filter.
d.PM <- d.PM[keep,]

############
## Obtain metric tensor from training data.
###########
porcent_train <- 0.5 #Proportion of data used to training.
samp_train <- c(seq(from = 1, to = round(dim(pickrell.eset)[2])*porcent_train), 
                seq(from = dim(pickrell.eset)[2] + 1, 
                 to = 1 + dim(pickrell.eset)[2] + round(dim(montgomery.eset)[2])*porcent_train)) 
d.PM_train <- DGEList(counts=d.PM$counts[,c(samp_train)], group = PM.eset$study[c(samp_train)])

group <- c(rep(1,round(dim(pickrell.eset)[2])*porcent_train), 
           rep(2,round(dim(montgomery.eset)[2])*porcent_train + 1))
counts.suvrel <-  Suvrel.normalization(d.PM_train$counts, group)
suvrel_results <- Suvrel.calc(counts.suvrel, group)

############
## MDS plot
###########
porcent_teste <- 1 - porcent_train  #Proportion of data used to test.
samp_teste <- 
  c(seq(from = round(dim(pickrell.eset)[2]*porcent_train) + 1, to = dim(pickrell.eset)[2]), 
    seq(from = 2 + dim(pickrell.eset)[2] + round(dim(montgomery.eset)[2]*porcent_train), 
        to = dim(pickrell.eset)[2] + round(dim(montgomery.eset)[2]))) 
d.PM_teste <- DGEList(counts=d.PM$counts[,c(samp_teste)], group = PM.eset$study[c(samp_teste)])

#Without feature selection
plot2file=TRUE
if(plot2file){
  setEPS()
  postscript(((paste("NoFestureSelection_", Sys.Date(),"_training_",porcent_train, "_.eps", sep=''))))
}

d.PM_teste_dat <- t(d.PM_teste$counts)
d <- dist(d.PM_teste_dat) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
fit # view results
x <- fit$points[,1]
y <- fit$points[,2]
k <- c(rep("green", dim(pickrell.eset)[2] - round(dim(pickrell.eset)[2]*porcent_train)) 
, rep("blue",(round(dim(pickrell.eset)[2]) + round(dim(montgomery.eset)[2])) 
      - (round(dim(pickrell.eset)[2]) + round(dim(montgomery.eset)[2]*porcent_train))))

########
#Mean and variance
########
pickrell_sample <- cbind(x[1:length(which(k=="green"))], y[1:length(which(k=="green"))])
pickrell_sample <- dist(pickrell_sample)
mean_pickrell_sample <- mean(pickrell_sample)
var_pickrell_sample <- sqrt(var(pickrell_sample))
message(paste("Mean, pick", mean_pickrell_sample, "Var, pick", var_pickrell_sample))

montgomery_sample <- cbind(x[(length(which(k=="green"))):(length(which(k=="green")) 
                                                        + length(which(k=="blue"))-1)], 
                         y[(length(which(k=="green"))):(length(which(k=="green")) 
                                                        + length(which(k=="blue"))-1)])
montgomery_sample <- dist(montgomery_sample)
mean_montgomery_sample <- mean(montgomery_sample)
var_montgomery_sample <- sqrt(var(montgomery_sample))
message(paste("Mean, montgomery", mean_montgomery_sample, "Var, montgomery", var_montgomery_sample))

plot(x, y, type="p", main = "MDS without feature selection", xlab="", ylab="", 
     lwd=10, col=k, pch=20, cex.main=1.85)
if(plot2file){
  dev.off()
}

#Feature selection
plot2file=TRUE
if(plot2file){
  setEPS()
  postscript(((paste("FestureSelection_", Sys.Date(),"_teste_",porcent_teste, "_.eps", sep=''))))
}

dat <- as.matrix(d.PM_teste$counts)
suv_sqrt <- sqrt(unlist(suvrel_results))
dat <- suv_sqrt*dat
dat <- t(dat)
d <- dist(dat) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
fit # view results
x <- fit$points[,1]
y <- fit$points[,2]
k <- c(rep("green", dim(pickrell.eset)[2] - round(dim(pickrell.eset)[2]*porcent_train)) 
       , rep("blue",(round(dim(pickrell.eset)[2]) + round(dim(montgomery.eset)[2])) 
             - (round(dim(pickrell.eset)[2]) + round(dim(montgomery.eset)[2]*porcent_train))))

########
#Mean and variance
########
pickrell_sample <- cbind(x[1:length(which(k=="green"))], y[1:length(which(k=="green"))])
pickrell_sample <- dist(pickrell_sample)
mean_pickrell_sample <- mean(pickrell_sample)
var_pickrell_sample <- sqrt(var(pickrell_sample))
message(paste("Mean, pick", mean_pickrell_sample, "Var, pick", var_pickrell_sample))

montgomery_sample <- cbind(x[(length(which(k=="green"))):(length(which(k=="green")) 
                                                          + length(which(k=="blue"))-1)], 
                           y[(length(which(k=="green"))):(length(which(k=="green")) 
                                                          + length(which(k=="blue"))-1)])
montgomery_sample <- dist(montgomery_sample)
mean_montgomery_sample <- mean(montgomery_sample)
var_montgomery_sample <- sqrt(var(montgomery_sample))
message(paste("Mean, montgomery", mean_montgomery_sample, "Var, montgomery", var_montgomery_sample))

plot(x, y, type="p", main = "MDS with feature selection", 
     xlab="", ylab="", lwd=10, col=k, pch=20, cex.main=1.85)

if(plot2file){
  dev.off()
}