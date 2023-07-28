##############
# libraries #
##############

library(adegenet)
library(poppr)
library(hierfstat)
library(PopGenReport)
library(pegas)
library(tidyverse)
############################
#     Load Data Files      #
############################
#load in genepop file as a genind object 
#change working directory 
setwd("C:/Users/laguiniga/Documents/GitHub/QUTO_Genetic_Analyses")

#load in data file
QUTO_genind <- read.genepop("QUTO/Data_Files/QUTO_byisland.gen",
                            ncode = 3)
