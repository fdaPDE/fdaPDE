# install stable version of fdaPDE

if("fdaPDE" %in% (.packages())){
  detach("package:fdaPDE", unload=TRUE) 
}

if("purrr" %in% (.packages())){
  detach("package:purrr", unload=TRUE) 
}

if("pracma" %in% (.packages())){
  detach("package:pracma", unload=TRUE) 
}

# create tests/tmp/ directory 
setwd("../")
dir.create("tests/tmp/")
install.packages("fdaPDE", lib = "tests/tmp/", repos="https://cran.mirror.garr.it/CRAN/")
install.packages("purrr", lib = "tests/tmp/")
install.packages("pracma", lib = "tests/tmp/")
library("fdaPDE", lib.loc = "tests/tmp/")
library("purrr", lib.loc = "tests/tmp/")
library("pracma", lib.loc = "tests/tmp/")

source("data/setup.R")

detach("package:fdaPDE", unload = TRUE)
remove.packages("fdaPDE", lib = "tests/tmp/")

# run tests
library(fdaPDE) # develop version
library(testthat)
test_dir("testthat/")

# delete tests/tmp directory
detach("package:purrr", unload=TRUE) 
detach("package:pracma", unload=TRUE)
remove.packages("purrr", lib = "tests/tmp/")
remove.packages("pracma", lib = "tests/tmp/")
unlink("tests/tmp/", recursive = TRUE)
setwd("tests/")