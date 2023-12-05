if(!require(pacman, quietly = TRUE)) install.packages("pacman")
pacman::p_load("pracma", "purrr", "testthat")

setwd("../")
# build data ------------------------------------------------------------------- 
if(!all(dir.exists(c(test_path("../data/DE-PDE"),
                     test_path("../data/SR-PDE"),
                     test_path("../data/tSR-PDE"),
                     test_path("../data/GSR-PDE"),
                     test_path("../data/tGSR-PDE"),
                     test_path("../data/fPCA"),
                     test_path("../data/STDE-PDE"))))){
  
  devtools::install_github("fdaPDE/fdaPDE", ref = "master")
  
  library(fdaPDE, quietly = TRUE)
  source("tests/data/setup.R")
  setwd("tests/")
  remove.packages("fdaPDE")
}