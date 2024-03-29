% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.FEM}
\alias{plot.FEM}
\title{Plot a \code{FEM} object}
\usage{
\method{plot}{FEM}(x, colormap = "heat.colors", num_refinements = NULL, ...)
}
\arguments{
\item{x}{A \code{FEM} object.}

\item{colormap}{A colormap exploited in the plot. The default value is the heat colormap.}

\item{num_refinements}{A natural number specifying how many bisections should be applied to each triangular element for
plotting purposes. This functionality is useful where a discretization with 2nd order Finite Element is applied.
This parameter can be specified only when a FEM object defined over a 2D mesh is plotted.}

\item{...}{Arguments representing graphical options to be passed to \link[rgl]{plot3d}.}
}
\value{
No return value
}
\description{
Three-dimensional plot of a \code{FEM} object, generated by \code{FEM} or returned by
\code{smooth.FEM} or \code{FPCA.FEM}.
If the \code{mesh} of the \code{FEMbasis} component is of class \code{mesh.2D} both the 3rd axis and the color represent
the value of the coefficients for the Finite Element basis expansion (\code{coeff} component of the \code{FEM} object).
If the \code{mesh} is of class \code{mesh.3D}, the color of each triangle or tetrahedron represent the mean value of
the coefficients for the Finite Element basis expansion (\code{coeff}).
}
\examples{
library(fdaPDE)
## Upload the horseshoe2D data
data(horseshoe2D)
boundary_nodes = horseshoe2D$boundary_nodes
boundary_segments = horseshoe2D$boundary_segments
locations = horseshoe2D$locations

## Create the 2D mesh
mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations), segments = boundary_segments)
## Create the FEM basis
FEMbasis = create.FEM.basis(mesh)
## Compute the coeff vector evaluating the desired function at the mesh nodes
## In this case we consider the fs.test() function introduced by Wood et al. 2008
coeff = fs.test(mesh$nodes[,1], mesh$nodes[,2])
## Create the FEM object
FEMfunction = FEM(coeff, FEMbasis)

## Plot the FEM function
plot(FEMfunction)
}
\seealso{
\code{\link{FEM}}, \code{\link{image.FEM}}
}
