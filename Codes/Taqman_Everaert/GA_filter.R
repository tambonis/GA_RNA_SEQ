################################################################################
################################################################################
# Objective: Filter the genomic feature that contains 0 counts.
# Author: Tiago Tambonis.
# Additional informations: 
# Date: 08/16. 
################################################################################
################################################################################

GA.filter <- function(counts.dat){
  filter <- as.logical(rowSums(counts.dat) > 0)
  counts.dat <- counts.dat[filter, ]
  return(counts.dat)
}
