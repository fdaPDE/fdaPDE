# fdaPDE

This repository contains the development version of fdaPDE package. In particular it contains the possibility of performing inference on the linear parameters, in exact and non exact way. Two possible inferential approaches are available for confidence intervals, three approaches for hypotesis testing.

##New features:
New features wrt CRAN: inference over linear parameters in exact way

FSPAI utility is now available for all the users

Compiled in Ubuntu using g++ compiler and in macOS: for the precise versions tested, see the report. If using a Linux machine, it is advisable to install rgl, geometry, plot3D and plot3Drgl before fdaPDE.

## Installation:
Two different methods are proposed in order to install the package in the R environment.  
Download the `.zip` file from the repository, unzip it, and for the installation choose one of the two following methods:  

- R console:
        ```install.packages("/path/to/fdaPDE-dev", type='source', repos=NULL)```

- From the Terminal: 
        ```$ R CMD build <path to folder to be installed>```     
        ```$ R CMD INSTALL -l <path name of the R library tree> <path to folder to be installed>```

see the installation section in the report for more information.
## Remarks:

1) the shift of indexes from R to C++ is done within the R functions smooth.FEM and FPCA.FEM Do not use C++ scripts directly on the R mesh objects, unless you take care of shifing indexes by yourself.
