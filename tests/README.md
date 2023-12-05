# Testing the development version of fdaPDE package

Compare the development version of the `fdaPDE` package with the stable version available in the `fdaPDE/fdaPDE` repository. 

### Install docker 
Docker is needed to compare the development version of `fdaPDE` package with the stable version available on CRAN.
To install docker, run from terminal:
```
  sudo apt-get install -y docker.io
  sudo usermod -aG docker $USER
```
Finally, reboot your PC or run from terminal:
```
  sudo reboot
```

### Run tests
Assuming your working directory is set to `fdaPDE/tests/` and you have installed 'docker' as explained in the previous section, to check the development version of the package, run the following command from the terminal:
```
  ./run_tests.sh
```
The bash script `run_tests.sh` will generate reference solutions, only once, using the stable version of the `fdaPDE` package. It will then compare these reference solutions with the outputs of your version of the package, exploiting the R package `testthat`. 
