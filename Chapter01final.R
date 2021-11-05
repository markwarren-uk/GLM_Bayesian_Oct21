
#=======================================

# R code for Chapter 1 of "Bayesian GLMs in R for Ecology" 
# by Mark Warren & Carl Smith

#=======================================

#Create an object 'a'
a <- 1
a

#Manipulate object
a + 1

#Another object
b <- a + 1
b

#Create vector 'freq'
freq <- c(40,99,31,35,0)
freq

class(freq)

#Create categorical vector 'species'
species <- c("pike", "roach", "chub", "perch", "asp")
species

class(species)

#Create a data frame from the vectors 'species' and 'freq'
spec_freq <- data.frame(species,freq) 

# Confirm the object contains the values by typing:
spec_freq

# What type of object is 'spec_freq'?
class(spec_freq)

#Import blue tit data

cyan <- read.table(file = "cyanistes.txt", 
                   header = TRUE, dec = ".")

#Install package "psych"
install.packages("psych")

#Load package "psych"
library(psych)

describe(cyan$depth, skew = FALSE)

#    vars   n mean  sd  min  max range se
# X1    1 438 0.33 0.1 0.17 0.75  0.58  0

#======================================= END
