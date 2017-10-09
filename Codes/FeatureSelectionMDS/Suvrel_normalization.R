################################################################################
##Mean 0, variance 1 normalization.
##Tiago Tambonis, 06/15.
################################################################################

Suvrel.normalization <- function(counts, group){      
    sd_g <- sqrt(apply(counts, 1, var))
    
    filter <- sd_g > 0 
    
    sd_g.min <- min(sd_g[filter])
    
    MEAN <- apply(counts, 1, mean)
    
    for (i in seq(dim(counts)[1]))
    {
      if (filter[i]==TRUE){
        counts[i,] <- (counts[i,] - MEAN[i])/sd_g[i]
      }else {(counts[i,]-MEAN[i])/sd_g.min}
    } 
  return(counts)
} 