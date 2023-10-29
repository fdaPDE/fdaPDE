# install stable version of fdaPDE
if(!require(pacman)) install.packages("pacman")
pacman::p_load("pracma", "purrr", "testthat")

if("fdaPDE" %in% (.packages())){
  detach("package:fdaPDE", unload=TRUE) 
}

# create tests/tmp/ directory 
setwd("../")
if(all(dir.exists(c(test_path("../data/DE-PDE"),
                    test_path("../data/SR-PDE"),
                    test_path("../data/tSR-PDE"),
                    test_path("../data/GSR-PDE"))))){
  dir.create("tests/tmp/")
  install.packages("fdaPDE", lib = "tests/tmp/", repos="https://cran.mirror.garr.it/CRAN/")

  library("fdaPDE", lib.loc = "tests/tmp/")
  source("tests/data/setup.R")

  detach("package:fdaPDE", unload = TRUE)
  remove.packages("fdaPDE", lib = "tests/tmp/")
  unlink("tests/tmp/", recursive = TRUE)
}

# run tests
library(fdaPDE) # develop version
test_dir("tests/testthat/")

# delete tests/tmp directory
detach("package:purrr", unload=TRUE) 
detach("package:pracma", unload=TRUE)
setwd("tests/")