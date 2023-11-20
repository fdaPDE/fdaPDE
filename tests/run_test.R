if(!require(pacman, quietly = TRUE)) install.packages("pacman")
pacman::p_load("pracma", "purrr", "testthat")

# install stable version of fdaPDE (build reference solutions only once) -------
if("fdaPDE" %in% (.packages())){
  detach("package:fdaPDE", unload=TRUE) 
}

# create tests/tmp/ directory 
setwd("../")
if(!all(dir.exists(c(test_path("../data/DE-PDE"),
                     test_path("../data/SR-PDE"),
                     test_path("../data/tSR-PDE"),
                     test_path("../data/GSR-PDE"),
                     test_path("../data/tGSR-PDE"),
                     test_path("../data/fPCA"),
                     test_path("../data/STDE-PDE")))){
  dir.create("tests/tmp/")
  install.packages("fdaPDE", lib = "tests/tmp/", repos="https://cran.mirror.garr.it/CRAN/")

  library("fdaPDE", lib.loc = "tests/tmp/", quietly = TRUE)
  
  source("tests/data/setup.R")

  detach("package:fdaPDE", unload = TRUE)
  remove.packages("fdaPDE", lib = "tests/tmp/")
  unlink("tests/tmp/", recursive = TRUE)
}

# run tests --------------------------------------------------------------------
library(fdaPDE, quietly = TRUE) # develop version
msg_file <- file("tests/msg.txt", open = "wt")
sink(msg_file, type = "message")

test_dir("tests/testthat/")

sink(type="message")
close(msg_file)
unlink("tests/msg.txt")

# delete tests/tmp directory
detach("package:purrr", unload=TRUE) 
detach("package:pracma", unload=TRUE)