# seventh test: Baumgartner's example "highdim" -- extremely large data set
library(cna)

setwd("..") # export output file to parent folder

sink(file = "r_output.txt") # start to export R output into file 

print(cna(d.highdim, # data set "highdim"
  rm.dup.factors=FALSE, # do not discard logically equivalent factors
  maxstep=c(5,5,10), # at most 5 conjuncts, 5 disjuncts and 10 factors per formula
  details = FALSE),
  nsolutions = 1) # return only one solution

sink(file = NULL) # stop exporting the output into file
