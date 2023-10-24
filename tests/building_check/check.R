install.packages(c('rgl', 'Matrix', 'plot3D', 'RcppEigen', 'Rcpp', 'MASS', 'testthat', 'rcmdcheck'), 
                 repos='https://cran.stat.unipd.it/')

rcmdcheck::rcmdcheck('fdaPDE/', check_dir= 'fdaPDE/fdaPDE.Rcheck/')