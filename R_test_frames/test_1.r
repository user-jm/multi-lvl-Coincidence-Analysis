# first test -- causal structure from Harbecke 2018 with an added upper level (A*B->C)
library(cna)

rnames <- c(1,2,3,4,5,6,7,8) # row names
cnames <- c("F1","F2","F3","F4","F5","F6","F7","A","B","C") # column names = names of causal factors

data_set <- array(0, dim=c(8,10),dimnames=list(rnames,cnames)) # a table of 8 rows and 10 columns each referenced by the entries in 
# rnames and cnames

# first row all factors are true
data_set[1,] <- 1

# second row F1,F3,F5,F6,B are true
true <- c("F1","F3","F5","F6","B")
data_set[2,true] <- 1

# third row F1,F2,F4,A are true
true <- c("F1","F2","F4","A")
data_set[3,true] <- 1

# fourth row only F1 is true
data_set[4,"F1"] <- 1

# fifth row F2,F3,F5,F6,B are true
true <- c("F2","F3","F5","F6","B")
data_set[5,true] <- 1

# sixth row F3,F5,F6,B are true
true <- c("F3","F5","F6","B")
data_set[6,true] <- 1

# seventh row only F2 is true
data_set[7,"F2"] <- 1

# eigth row all factors are false

setwd("..") # export output file to parent folder

sink(file = "r_output.txt") # start to export R output into file 

print(cna(data_set, # use the truth table data_set
  rm.dup.factors=FALSE, # do not discard logically equivalent factors
  maxstep=c(5,5,10), # at most 5 conjuncts, 5 disjuncts and 10 factors per formula
  details = FALSE,
  ordering = list(c("F1","F2","F3","F4","F5","F6","F7"),c("A","B","C"))), # ordering(e_1,e_2,...) places e_2 causally downstream to e_1
  # this is how the constitution levels get separated
  nsolutions = "all") # return all solutions

sink(file = NULL) # stop exporting into file
