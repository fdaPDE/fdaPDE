# Testing the development version of fdaPDE package

Compare the development version of `fdaPDE` package with the stable version available on CRAN.  
To run correctly the tests, check whether the following conditions holds:

  - your working directory is set to `fdaPDE/tests/`.
  
  - your development version of `fdaPDE` is loaded through `library(fdaPDE)`
  
Finally, to check the development version of the package, run from terminal:
```
  Rscript run_test.R
```
The script `run_test.R` will generate, only once, reference solutions relying on the stable version of `fdaPDE` package and it will compare these reference solutions with the outputs of your version of the package exploiting the R package `testthat`. 
