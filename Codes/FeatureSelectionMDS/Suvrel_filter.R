################################################################################
##This function executes the filtering of the counts.
##The genes that will go away are that sum of column are less or equal to zero.
##Tiago Tambonis, 06/14.
################################################################################

Suvrel.filter <- function(counts){
  filter <- as.logical(rowSums(counts) > 0)
  counts <- counts[filter, ]
  return(counts)
}
