test_that("tGSR-PDE 2d",{
  
  options(warn=-1)
  foldername <- test_path("../data/tGSR-PDE/test_25/")
  
  FAMILY = "gamma"
  inv.link <-  function(mu){ return(-1 / mu)} 
  
  data("horseshoe2D")
  mesh = create.mesh.2D(nodes = horseshoe2D$boundary_nodes,
                        segments = horseshoe2D$boundary_segments)
  mesh = refine.mesh.2D(mesh, maximum_area =  0.025, minimum_angle = 30 )
  FEMbasis = create.FEM.basis(mesh)
  
  f <- function(x, y, t) {
    
    return((-1.5 * sin(2 * pi * x) * cos(2 * pi * y) +
              2 / 5 * sin(3 * pi * x * t) - 2) * (t + 1))
  }
  
  is.in.horseshoe <- function(x, y) {
    r <- .5
    r0 <- .1
    
    if ((x - 3)^2 + (y - r)^2 < (r - r0)^2 & x > 3) {
      return(TRUE)
    }
    if (r0^2 <= x^2 + y^2 & x^2 + y^2 <= (2 * r - r0)^2 & x <= 0) {
      return(TRUE)
    }
    if ((x - 3)^2 + (y + r)^2 < (r - r0)^2 & x > 3) {
      return(TRUE)
    }
    if (abs(y) > r0 & abs(y) < 2 * r - r0 & 0 < x & x <= 3) {
      return(TRUE)
    }
    return(FALSE)
  }
  
  is.p.in.horseshoe <- function(p) {
    return(is.in.horseshoe(p[1], p[2]))
  }
  
  ### generating data ###
  M = 4
  m = M
  n = 400
  
  time_mesh <- seq(0, 1, length.out = M)
  time_locations = time_mesh
  
  set.seed(32)
  
  loc = cbind(runif(2 * n, min = -1, max = 4),
              runif(2 * n, min = -1, max = 1))
  ww <- apply(loc, 1, is.p.in.horseshoe) # ! is.na(fs.test(loc[,1],
  # loc[ ,2], exclude = T))
  loc <- loc[ww, ]  
  
  space_time_locations <-cbind(rep(time_locations, each = nrow(loc)),
                               rep(loc[, 1], length(time_locations)),
                               rep(loc[, 2], length(time_locations)) )
  ### covariates ###
  betas = c(-.2, .3)
  desmat <- matrix(0, nrow = nrow(space_time_locations), ncol = 2)
  desmat[, 1] <- rbeta(n = nrow(desmat), shape1 = 1.5, shape2 = 2)
  desmat[, 2] <- rbeta(n = nrow(desmat), shape1 = 3, shape2 = 2)
  
  true.field <- f(space_time_locations[, 2], space_time_locations[, 3],space_time_locations[, 1])  
  field = true.field + desmat %*% betas 
  
  scale = 1
  data <- rgamma(n = nrow(space_time_locations),
                 shape = -1 / field / scale, scale = scale)
  data <- matrix(data, nrow = nrow(loc), ncol = m)
  
  ### smoothing parameters ###
  lambdaS <- 10^seq(-3, -1, by=0.5)
  lambdaT <- 10^seq(-3, -2, by=0.5)
  
  ########### PARABOLIC ###########
  
  #### Test 1.1: without GCV
  output_CPP <- smooth.FEM.time(time_mesh = time_mesh, time_locations = time_locations,
                                locations = loc, observations = data,  FEMbasis = FEMbasis, 
                                covariates = desmat,  
                                lambdaS = lambdaS[2], lambdaT = lambdaT[2],
                                family = FAMILY, FLAG_PARABOLIC = T)
  
  #### Test 1.2: with stochastic GCV
  output_CPP <- smooth.FEM.time(time_mesh = time_mesh, time_locations = time_locations,
                                locations = loc, observations = data,  FEMbasis = FEMbasis, 
                                covariates = desmat,  
                                DOF.evaluation = "stochastic", 
                                lambda.selection.lossfunction = "GCV",
                                lambdaS = lambdaS, lambdaT = lambdaT,
                                family = FAMILY, FLAG_PARABOLIC = T)
  
  ########### SEPARABLE ############
  
  #### Test 1.3: without GCV
  output_CPP <- smooth.FEM.time(time_mesh = time_mesh, time_locations = time_locations,
                                locations = loc, observations = data, FEMbasis = FEMbasis, 
                                covariates = desmat, 
                                lambdaS = lambdaS[2], lambdaT = lambdaT[2],
                                family = FAMILY, FLAG_PARABOLIC = F)
  
  #### Test 1.4: with stochastic GCV
  output_CPP <- smooth.FEM.time(time_mesh = time_mesh, time_locations = time_locations,
                                locations = loc, observations = data, FEMbasis = FEMbasis, 
                                covariates = desmat, 
                                DOF.evaluation = "stochastic", lambda.selection.lossfunction = "GCV",
                                lambdaS = lambdaS, lambdaT = lambdaT,
                                family = FAMILY, FLAG_PARABOLIC = F)
  
  
  FAMILY = "gamma"
  inv.link <-  function(mu){ return(-1 / mu)} 
  
  data("horseshoe2D")
  mesh = create.mesh.2D(nodes = horseshoe2D$boundary_nodes,
                        segments = horseshoe2D$boundary_segments)
  mesh = refine.mesh.2D(mesh, maximum_area =  0.025, minimum_angle = 30 )
  FEMbasis = create.FEM.basis(mesh)
  
  f <- function(x, y, t) {
    
    return((-1.5 * sin(2 * pi * x) * cos(2 * pi * y) +
              2 / 5 * sin(3 * pi * x * t) - 2) * (t + 1))
  }
  
  is.in.horseshoe <- function(x, y) {
    r <- .5
    r0 <- .1
    
    if ((x - 3)^2 + (y - r)^2 < (r - r0)^2 & x > 3) {
      return(TRUE)
    }
    if (r0^2 <= x^2 + y^2 & x^2 + y^2 <= (2 * r - r0)^2 & x <= 0) {
      return(TRUE)
    }
    if ((x - 3)^2 + (y + r)^2 < (r - r0)^2 & x > 3) {
      return(TRUE)
    }
    if (abs(y) > r0 & abs(y) < 2 * r - r0 & 0 < x & x <= 3) {
      return(TRUE)
    }
    return(FALSE)
  }
  
  is.p.in.horseshoe <- function(p) {
    return(is.in.horseshoe(p[1], p[2]))
  }
  
  ### generating data ###
  M = 4
  m = M
  n = 400
  
  time_mesh <- seq(0, 1, length.out = M)
  time_locations = time_mesh
  
  set.seed(32)
  
  loc = cbind(runif(2 * n, min = -1, max = 4),
              runif(2 * n, min = -1, max = 1))
  ww <- apply(loc, 1, is.p.in.horseshoe) # ! is.na(fs.test(loc[,1],
  # loc[ ,2], exclude = T))
  loc <- loc[ww, ]  
  
  space_time_locations <-cbind(rep(time_locations, each = nrow(loc)),
                               rep(loc[, 1], length(time_locations)),
                               rep(loc[, 2], length(time_locations)) )
  ### covariates ###
  betas = c(-.2, .3)
  desmat <- matrix(0, nrow = nrow(space_time_locations), ncol = 2)
  desmat[, 1] <- rbeta(n = nrow(desmat), shape1 = 1.5, shape2 = 2)
  desmat[, 2] <- rbeta(n = nrow(desmat), shape1 = 3, shape2 = 2)
  
  true.field <- f(space_time_locations[, 2], space_time_locations[, 3],space_time_locations[, 1])  
  field = true.field + desmat %*% betas 
  
  scale = 1
  data <- rgamma(n = nrow(space_time_locations),
                 shape = -1 / field / scale, scale = scale)
  data <- matrix(data, nrow = nrow(loc), ncol = m)
  
  ### smoothing parameters ###
  lambdaS <- 10^seq(-3, -1, by=0.5)
  lambdaT <- 10^seq(-3, -2, by=0.5)
  
  #### Test 25.1: without GCV (PARABOLIC)
  invisible(capture.output(sol <- smooth.FEM.time(time_mesh = time_mesh, time_locations = time_locations,
                                                      locations = loc, observations = data,  FEMbasis = FEMbasis, 
                                                      covariates = desmat,  
                                                      lambdaS = lambdaS[2], lambdaT = lambdaT[2],
                                                      family = FAMILY, FLAG_PARABOLIC = T)))
  load(file=paste0(foldername,"/test_25_1.RData"))
  expect_equal( max(abs((sol$fit.FEM.time$coeff-sol_ref$fit.FEM.time$coeff))) < tol, TRUE);
  
  #### Test 25.2: without GCV (SEPARABLE)
  invisible(capture.output(sol <- smooth.FEM.time(time_mesh = time_mesh, time_locations = time_locations,
                                                      locations = loc, observations = data, FEMbasis = FEMbasis, 
                                                      covariates = desmat, 
                                                      lambdaS = lambdaS[2], lambdaT = lambdaT[2],
                                                      family = FAMILY, FLAG_PARABOLIC = F)))
  load(file=paste0(foldername,"/test_25_2.RData"))
  expect_equal( max(abs((sol$fit.FEM.time$coeff-sol_ref$fit.FEM.time$coeff))) < tol, TRUE);
})