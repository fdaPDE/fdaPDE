if(!require(pacman, quietly = TRUE)) install.packages("pacman")
pacman::p_load("pracma", "purrr", "testthat")

# run tests --------------------------------------------------------------------
setwd("../")

# install develop version
install.packages(".", repos = NULL, type = "source") 
# load develop version
library(fdaPDE, quietly = TRUE)

msg_file <- file("tests/msg.txt", open = "wt")
sink(msg_file, type = "message")

test_dir("tests/testthat/")

sink(type="message")
close(msg_file)
unlink("tests/msg.txt")
