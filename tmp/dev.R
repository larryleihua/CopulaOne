# for developing the CopulaOne package
rm(list = setdiff(ls(), lsf.str()))
library(devtools)
library(roxygen2)
setwd("/home/lh/Dropbox/git/CopulaOne")
#setwd("E:/Dropbox/git/CopulaOne")

document()
setwd("..")
install("CopulaOne")
library(CopulaOne)

# exclude folders for building bundle
# devtools::use_build_ignore("notes")
# library(rticles) # rmarkdown template

### add forex data ###
USDCAD2015_H1 = read.csv("~/Dropbox/justdoit/TD049-TailMomentum/R/USDCAD_2015_H1_NEWYORK.csv", sep=",")
AUDUSD2015_H1 = read.csv("~/Dropbox/justdoit/TD049-TailMomentum/R/AUDUSD_2015_H1_NEWYORK.csv", sep=",")
use_data(USDCAD2015_H1)
use_data(AUDUSD2015_H1)

