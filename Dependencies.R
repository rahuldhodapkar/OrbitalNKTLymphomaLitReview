#!/usr/bin/env Rscript
## Dependencies.R
##
# Small utility script to consolidate the requirements for data analysis
# and ensure easy installation for replication and extension of findings.
#
# Note that on macOS, additional dependencies are required, including
# - Xcode command line utilities
# - `brew install openssl curl libgit2`
#
# @author Rahul Dhodapkar <rahul.dhodapkar@yale.edu>
# @version 2021.02.10
#
# Last tested on 2022.01.28
#

## set local CRAN mirror
local({r <- getOption("repos")
       r["CRAN"] <- "http://cran.r-project.org" 
       options(repos=r)
})

##############################################
# INSTALL FROM CRAN
##############################################

install.packages("ggplot2")
install.packages("survival")
install.packages("survminer")
install.packages("stringr")
install.packages("cowplot")
install.packages("tidyverse")
install.packages("caret")
install.packages("leaps")
