# eighth test -- a frame to generate CNA-output directly from a formula
# to be adjusted:
# formula = the generating formula
# the values of maxstep=c(5,5,10) in cna(...) (see comments below)

# syntax:
# negation: minuscle of factor (not(A)=a)
# conjunction: "*"
# disjunction: "+"
# implication: "->"
# equivalence: "<->"
library(cna)


formula <-"(A*B <-> C)*(C*A*d <-> E)"  # condition to generate the truth table

data_set <-selectCases(formula) # create truth table


print(data_set)  # print table into console

setwd("..") # export output file to parent folder

sink(file = "r_output.txt") # start to export R output into file 

print(cna(data_set, # use the truth table data_set
  rm.dup.factors=FALSE, # do not discard logically equivalent factors
  maxstep=c(5,5,10), # at most 5 conjuncts, 5 disjuncts and 10 factors per formula
  details = FALSE,
  # for multi-level structures: comment out the following line and modify it accordingly
  #ordering = list(c("V1","V2","V3","V4","V5","V6","V7"),c("V8","V9","V10")) # ordering(e_1,e_2,...) places e_2 causally downstream to e_1
  # this is how the constitution levels get separated
  ), nsolutions = "all") # return all solutions

sink(file = NULL) # stop exporting into file
