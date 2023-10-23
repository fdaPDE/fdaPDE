dir.create("tests/tmp")
install.packages("fdaPDE", lib = "tests/tmp/", repos="https://cran.mirror.garr.it/CRAN/")
library("fdaPDE", lib.loc = "tests/tmp")
library(testthat)
source("tests/data/setup.R")

detach("package:fdaPDE", unload = TRUE)
remove.packages("fdaPDE", lib = "tests/tmp/")
unlink("tests/tmp", recursive = TRUE)

library(fdaPDE) # develop version
library(testthat)
test_dir("tests/testthat/")
