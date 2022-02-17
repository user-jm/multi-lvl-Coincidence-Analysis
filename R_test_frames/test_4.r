# fourth test: relation with negations, simple single level structure
library(cna)

rnames <- c(1,2,3,4,5,6,7,8) # row names
cnames <- c("A","B","C","I") # column names = names of the causal factors

data_set <- array(0, dim=c(8,4),dimnames=list(rnames,cnames))  # a table of 8 rows and 4 columns each referenced by the entries in 
# rnames and cnames

# A is a necessary condition for I
true <- c("A","I")
data_set[1,true] <- 1

# necessary condition is satisfied, veto only partially
true <- c("A","B","I")
data_set[2,true] <- 1

# necessary condition is satisfied, veto only partially
true <- c("A","C","I")
data_set[3,true] <- 1

# A is true but, B*C vetoes
true <- c("A","B","C")
data_set[4,true] <- 1

# only B, C are true
true <- c("B","C")
data_set[5,true] <- 1

# only B is true
true <- c("B")
data_set[6,true] <- 1

# only C is true
true <- c("C")
data_set[7,true] <- 1

# nothing is true -> all values are zero

setwd("..") # export output file to parent folder

sink(file = "r_output.txt") # start to export R output into file 

print(cna(data_set, # use the truth table data_set
  rm.dup.factors=FALSE, # do not discard logically equivalent factors
  maxstep=c(3,3,4), # at most 3 conjuncts, 3 disjuncts and 4 factors per formula
  details = FALSE),
  nsolutions = "all") # return all solutions

sink(file = NULL) # stop exporting into file
