################################################################################
################################################################################
# Objective: Geometric approach. 
# Author: Tiago Tambonis and Marcelo Boareto.
# Additional informations: 
# Date: 08/2017. 
################################################################################
################################################################################

GA<- function(counts.dat, group){
  colnames(counts.dat) <- factor(group)

  condA <- counts.dat[,colnames(counts.dat)==1]
  condB <- counts.dat[,colnames(counts.dat)==2]
  
  f_diff <- data.frame()
  f_diff_temp = 0
  n_diff = 0
  for (k in seq(1, dim(counts.dat)[1])){
    for (i in seq(1, dim(condA)[2])){
      for (j in seq(1, dim(condB)[2])){
        f_diff_temp = f_diff_temp + abs(condA[k,i] - condB[k,j])
        n_diff = n_diff + 1
      }
    }
    f_diff = rbind(f_diff, f_diff_temp)
    f_diff_temp = 0  
  }
  
  n_same = 0
  f_same_A = data.frame()
  f_same_temp = 0 
  combinations <-  combn(1:(dim(condA)[2]), 2)
  for (k in seq(1, dim(counts.dat)[1])){
    for (i in seq(1:dim(combinations)[2])){
      f_same_temp = f_same_temp + abs(condA[k,combinations[1,i]] - condA[k,combinations[2,i]]) 
      n_same = n_same + 1 
    }
    f_same_A = rbind(f_same_A, f_same_temp)
    f_same_temp = 0 
  }

  f_same_B = data.frame()  
  f_same_temp = 0 
  combinations <-  combn(1:(dim(condB)[2]), 2)
  for (k in seq(1, dim(counts.dat)[1])){
    for (i in seq(1:dim(combinations)[2])){
      f_same_temp = f_same_temp + abs(condB[k,combinations[1,i]] - condB[k,combinations[2,i]]) 
      n_same = n_same + 1 
    }
    f_same_B = rbind(f_same_B, f_same_temp)
    f_same_temp = 0 
  }
  
  f_same = f_same_A + f_same_B
    
  w = (f_diff/n_diff)**2 - (f_same/n_same)**2
  
  w[w<0] = 0
  
  w = w/sqrt(sum((w**2)))
  
  rownames(w) <- rownames(counts.dat)
  colnames(w) <- "Results"
  
  return(w)
}
