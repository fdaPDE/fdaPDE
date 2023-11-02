# fdaPDE

This repository contains the development version of fdaPDE package. 

New features wrt CRAN: Wald inference on the nonparametric term for space-time models.

Compiled in Ubuntu using g++ compiler and in macOS: for the precise versions tested, see the report. If using a Linux machine, it is advisable to install rgl, geometry, plot3D and plot3Drgl before fdaPDE.

Image.FEM has been restored. Bugs in fPCA, boundary conditions and space-varying regression have been fixed. Issues of point location in 2.5D have been fixed.
Compiled in Win RStudio, Ubuntu using g++ compiler and in macOS: for the precise versions tested, see the report. If using a Linux machine, it is advisable to install rgl, geometry, plot3D and before fdaPDE. If using Windows, it is advisable to install Rtools and then rgl, plot3D, geometry and RcppEigen libraries. 

## Subfolder structure:
/src contains all C++ code and a special file named Makevars necessary to build and install the R package. The code is now organized in subfolders and divided into source files and header files. See the report for the precise new orgaization. We suggest to use a base-8 scale as spacing scale for visualization, in order to preserve the code style.

/R contains the R functions that wrap the C++ calls,

/data contains the data to run the tests in the examples,

/man contains the manual for the package.

## Installation:
Two different methods are proposed in order to install the package in the R environment.  
Download the `.zip` file from the repository, unzip it, and for the installation choose one of the two following methods:  

- R console:
        ```install.packages("/path/to/fdaPDE", type='source', repos=NULL)```

- From the Terminal: 
        ```$ R CMD build <path to folder to be installed>```     
        ```$ R CMD INSTALL -l <path name of the R library tree> <path to folder to be installed>```

see the installation section in the report for more information.

## Remarks:
1) the shift of indexes from R to C++ is done within the R functions smooth.FEM and FPCA.FEM Do not use C++ scripts directly on the R mesh objects, unless you take care of shifing indexes by yourself.
