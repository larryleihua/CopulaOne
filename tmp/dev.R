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
