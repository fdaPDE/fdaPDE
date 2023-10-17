
test_that("SR-PDE 2.5D Horseshoe domain",{
  
options(warn=-1)
foldername = test_path("../data/SR-PDE/test_11/")

data(horseshoe2.5D)
mesh = horseshoe2.5D

nnodes=dim(mesh$nodes)[1]

FEMbasis=create.FEM.basis(mesh)

# Exact solution (pointwise at nodes)
sol_exact=fs.test.3D(mesh$nodes[,1],mesh$nodes[,3],mesh$nodes[,2])

# Add error to simulate data
set.seed(7893475)
ran=range(sol_exact)
data = sol_exact + rnorm(nnodes, mean=0, sd=0.05*abs(ran[2]-ran[1]))

# Set smoothing parameter
lambda = 10^seq(-2,0.5,by=0.25)

#### Test 1.1: Without GCV
invisible(capture.output(sol<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda)))

load(file=paste0(foldername,"/test_11_1.RData"))
expect_equal( max(abs((sol$fit.FEM$coeff-output_CPP$fit.FEM$coeff))) < tol, TRUE);
#### Test 1.2: grid with exact GCV
#it takes a lot of time
invisible(capture.output(sol<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')))

load(file=paste0(foldername,"/test_11_2.RData"))
expect_equal( max(abs((sol$fit.FEM$coeff-output_CPP$fit.FEM$coeff))) < tol, TRUE);
})