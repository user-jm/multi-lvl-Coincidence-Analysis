# second test -- a simple causal chain on two levels
library(cna)

rnames <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16) # row names
cnames <- c("A","B","C","D","E","F","G","H","I","J") # column names = names of causal factors

data_set <- array(0, dim=c(16,10),dimnames=list(rnames,cnames)) # a table of 16 rows and 10 columns each referenced by the entries in 
# rnames and cnames

# first row all factors are true
data_set[1,] <- 1

# second row A,D,E,F,G are true
true <- c("A","D","E","F","G")
data_set[2,true] <- 1

# third row D,E,H are true
true <- c("A","D","E","F","H")
data_set[3,true] <- 1

# fourth row D,G,H are true
true <- c("B","D","G","H","I")
data_set[4,true] <- 1

# fifth row E,G,H are  true
true <- c("B","E","G","H","I")
data_set[5,true] <- 1

# sixth row D,E are  true
true <- c("A","D","E","F")
data_set[6,true] <- 1

# seventh row D,G are  true
true <- c("D","G")
data_set[7,true] <- 1

# eighth row D,H are  true
true <- c("D","H")
data_set[8,true] <- 1

# ninth row E,G are  true
true <- c("E","G")
data_set[9,true] <- 1

# tenth row E,H are  true
true <- c("E","H")
data_set[10,true] <- 1

# eleventh row H,G are  true
true <- c("B","H","G","I")
data_set[11,true] <- 1

# twelfth row D is  true
true <- c("D")
data_set[12,true] <- 1

# thirdteenth row E is  true
true <- c("E")
data_set[13,true] <- 1

# fourteenth row G is  true
true <- c("G")
data_set[14,true] <- 1

# fifteenth row H is  true
true <- c("H")
data_set[15,true] <- 1

# sixteenth row all factors are false

setwd("..") # export output file to parent folder

sink(file = "r_output.txt") # start to export R output into file 

print(cna(data_set, # use the truth table data_set
  rm.dup.factors=FALSE, # do not discard logically equivalent factors
  maxstep=c(5,5,10), # at most 5 conjuncts, 5 disjuncts and 10 factors per formula
  details = FALSE,
  ordering = list(c("D","E","F","G","H","I","J"),c("A", "B","C"))), # ordering(e_1,e_2,...) places e_2 causally downstream to e_1
  # this is how the constitution levels get separated
  nsolutions = "all") # return all solutions

sink(file = NULL) # stop exporting into file
