################################################################################
################################################################################
# Objective: RPM normalization.
# Author: Tiago Tambonis.
# Additional informations: 
# Date: 08/17. 
################################################################################
################################################################################

RPM_normalization <- function(counts.dat){
  ##################
  #Normalization
  ##################
  counts.dat <- counts.dat + 1 #Avoid problems with log scale.
  pm_factor <- apply(counts.dat, 2, function(x) sum(x)/1000000)
  counts.dat <- sweep(counts.dat,MARGIN=2,pm_factor,'/')
  
  return(counts.dat)
}
