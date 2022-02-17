# eighth test -- a frame to generate CNA-output directly from a formula
# to be adjusted:
# formula = the generating formula
# complete_table = the number of causal factors
# names(complete_table) = the names of the causal factors
# 
# the values of maxstep=c(5,5,10) in cna(...) (see comments below)

# syntax:
# negation: minuscle of factor (not(A)=a)
# conjunction: "*"
# disjunction: "+"
# implication: "->"
# equivalence: "<->"

library(cna)
library('data.table')

formula <-"(A*B <-> C)*(C*A*d <-> E)*(E <-> F)*(A*B*d <-> F)"  # condition to generate the truth table

complete_table <- allCombs(c(2,2,2,2,2,2)) - 1 # creates a table with all combinations of values 0 and 1 in each column
# (first the table contains the combinations of 1 and 2, but -1 reduces all entries
# make sure that there is one column for each variable 

names(complete_table) <- c("A","B","C","D","E","F") # specify the column names = names of causal factors

data_set <- selectCases(formula, complete_table)  # generate the truth table for formula

print(data_set)  # print the table into the console

setwd("..") # export output file to parent folder

sink(file = "r_output.txt") # start to export R output into file 

print(cna(data_set, # use the truth table data_set
  rm.dup.factors=FALSE, # do not discard logically equivalent factors
  maxstep=c(5,5,10), # at most 5 conjuncts, 5 disjuncts and 10 factors per formula
  details = FALSE,
  # for multi-level structures:
  ordering = list(c("A","B","C","D","E"),c("F")) # ordering(e_1,e_2,...) places e_2 causally downstream to e_1
  # this is how the constitution levels get separated
  ), nsolutions = "all") # return all solutions

sink(file = NULL) # stop exporting into file
