# Jens' Beispiel
library(cna)

rnames <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
cnames <- c("A","B","C","D","E","F","G","H","I","J")

data_set <- array(0, dim=c(16,10),dimnames=list(rnames,cnames))

# erste Zeile alles auf 1
data_set[1,] <- 1

# zweite Zeile D,E,G auf 1
wahr <- c("A","D","E","F","G")
data_set[2,wahr] <- 1

# dritte Zeile D,E,H auf 1
wahr <- c("A","D","E","F","H")
data_set[3,wahr] <- 1

# vierte Zeile D,G,H auf 1
wahr <- c("B","D","G","H","I")
data_set[4,wahr] <- 1

# fuenfte Zeile E,G,H auf 1
wahr <- c("B","E","G","H","I")
data_set[5,wahr] <- 1

# sechste Zeile D,E auf 1
wahr <- c("A","D","E","F")
data_set[6,wahr] <- 1

# siebente Zeile D,G auf 1
wahr <- c("D","G")
data_set[7,wahr] <- 1

# achte Zeile D,H auf 1
wahr <- c("D","H")
data_set[8,wahr] <- 1

# neunte Zeile E,G auf 1
wahr <- c("E","G")
data_set[9,wahr] <- 1

# zehnte Zeile E,H auf 1
wahr <- c("E","H")
data_set[10,wahr] <- 1

# elfte Zeile H,G auf 1
wahr <- c("B","H","G","I")
data_set[11,wahr] <- 1

# zwoelfte Zeile D auf 1
wahr <- c("D")
data_set[12,wahr] <- 1

# dreizehnte Zeile E auf 1
wahr <- c("E")
data_set[13,wahr] <- 1

# vierzehnte Zeile G auf 1
wahr <- c("G")
data_set[14,wahr] <- 1

# fuenfzehnte Zeile H auf 1
wahr <- c("H")
data_set[15,wahr] <- 1

# sechzehnte Zeile alles auf 0

setwd("..") # Ausgabe soll in uebergeordneten Ordner erzeugt werden

sink(file = "r_output.txt") # Ausgabe werden ab hier in Datei gespeichert

print(cna(data_set, # unserer Datensatz
  rm.dup.factors=FALSE, # verwerfe Spalten mit identischen Eintraegen nicht
  maxstep=c(5,5,10), # maximal 5 Konjunkte, 5 Disjunkte und 10 Faktoren
#  what = "a", # zeige nur atomare Loesungsformeln
  details = FALSE,
  strict = FALSE,
  ordering = list(c("D","E","F","G","H","I","J"),c("A", "B","C"))), # ordering(e_1,e_2,...) setzt e_2 downstream bezueglich e_1, hier werden die Ebenen separiert
  nsolutions = "all") # gib alle Loesungen an

sink(file = NULL) # Ausgabe in Datei endet hier
