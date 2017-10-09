################################################################################
##Suvrel implementation.
##Tiago Tambonis, 02/16.
################################################################################

Suvrel.calc <- function(counts, group){   
    
  colnames(counts) <- factor(group) #Imposition due the structure of calc_mean function.
                                        
  calc_mean <- function(counts, condition){ #Function for average calculation.
    mean <- rowMeans(counts[,group==condition])
    return(mean)
  }
  
  means <- sapply(seq(2), function(i) calc_mean(counts,i)) #Average calculation.

  epsilon.k <- apply(means, 1, function(x) -2.0*((x[1] - x[2])^2))
  
  relevances <- (-epsilon.k)/(sqrt(sum(epsilon.k^2))) #Calculation of relevances.
  
  results <- list(Relevances=relevances) #Output the relevances.
    
  return(results)
}
