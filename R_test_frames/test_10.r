# fifth test: Baumgartner's example "educate", simple single level structure
library(cna)

setwd("..") # export output file to parent folder

sink(file = "r_output.txt") # start to export R output into file 

print(cna(d.educate, # data set "educate"
  rm.dup.factors=FALSE, # do not discard logically equivalent factors
  maxstep=c(5,5,10), # at most 5 conjuncts, 5 disjuncts and 10 factors per formula
  details = FALSE),
  nsolutions = "all") # return all solutions

sink(file = NULL) # stop exporting the output into file
