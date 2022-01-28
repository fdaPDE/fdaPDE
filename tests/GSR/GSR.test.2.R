# Geometry + Spatial Field
{
  library(purrr)
  source("~/Scrivania/fdaPDE/tests/Auxiliary/Regressione_Mattina.R")
  source("~/Scrivania/fdaPDE/tests/Auxiliary/integrate_f.R")
  source("~/Scrivania/fdaPDE/tests/Auxiliary/isInside.R")
  source("tests/Auxiliary/R_plot_graph.ggplot2.R")
  source("~/Scrivania/fdaPDE/tests/Auxiliary/Eval_split_discontinous.R")
  source("~/Scrivania/fdaPDE/tests/Auxiliary/Dijkstra.R")
  
  #FAMILY = "binomial"
 # logit <- function(x){qlogis(x)}
#  inv.logit <- function(x){plogis(x)}
  #link = function(x) { log(x/(1-x))}
  #inv.link = function(x){ exp(x)/(1+exp(x))}
  FAMILY = "poisson"
  l<-make.link("log")
  link<-l$linkfun
  inv.link<-l$linkinv
  
  data("spiders")
  vertices = cbind(spiders$domain$vertices$x, spiders$domain$vertices$y)
  edges = cbind(spiders$domain$from, spiders$domain$to)
  
  L = as.linnet(spiders$domain)
  
  #580 nodes -> delta = 40 , delta = 50 -> 408
  delta =40
  mesh = create.mesh.1.5D(vertices, edges)
  mesh=refine.mesh.1.5D(mesh, delta=50)
  nnodes = nrow(mesh$nodes)
  FEMbasis = create.FEM.basis(mesh)
  
  plot(mesh)
  points(vertices[c(55,98),], pch=16)
  dijkstra.55 = Dijkstra(mesh,55)
  #dijkstra.66= Dijkstra(mesh,66)
  dijkstra.98 = Dijkstra(mesh, 98)
  
  AUX = function(x,y,seg,tp){
    sigma.55 = 75
    res.55 = equal_split_discontinous(mesh, sigma.55, dijkstra.55, x, y)
    
    sigma.98 = 75
    res.98= equal_split_discontinous(mesh, sigma.98, dijkstra.98, x, y)
    
    res = res.55$coef + res.98$coef 
    
    # idx.55 = which(res.55$bandwidth==0.0)
    # idx.98 = which(res.98$bandwidth==0.0)
    # 
    # intersect.tot = intersect(idx.55, idx.98)
    # res[intersect.tot] = -5.0
    
    idx.55 = which(res.55$bandwidth==1.0)
    idx.98 = which(res.98$bandwidth==1.0)
    
    res[idx.55] = 1000*res[idx.55] 
    res[idx.98] = 1000*res[idx.98] 
    
    
    return(res)
  }
  
  coef.ex = AUX(mesh$nodes[,1], mesh$nodes[,2])
  R_plot_graph.a.sym.ggplot2(FEM(coef.ex, FEMbasis))
}

{
  param = coef.ex
  mu<-inv.link(param)
  range(mu)
  # sampling response:
  response <- rpois(nrow(mesh$nodes), lambda = mu)
  
  # palette = c("green4", "red3")
  # colors = c()
  # for(i in 1:nrow(mesh$nodes)){
  #   if( response[i] == TRUE )
  #     colors = append(colors, palette[1])
  #   else
  #     colors = append(colors, palette[2])
  # }
  heat <- heat.colors(n=24)
  p = heat[round(99*(response-min(response))/diff(range(response)))+1]
  
  plot(mesh)
  points(mesh$nodes, pch=16, col=p)
  
}


nobs = c(200, 400, 600, 800, 1000, 1200)
nobs=c(600, 800)

N = length(nobs)
M = 1
times = matrix(0, nrow=M, ncol=N)
mise = matrix(0,nrow=M, ncol=N)
rmse = matrix(0, nrow=M, ncol=N)
lambda = 10^seq(1,2,length.out = 20)
beta_ex =1.0

for(i in 1:N){
  for(j in 1:M){
    PP = runiflpp(nobs[i], L)
    locations = cbind(PP$data$x, PP$data$y)
    
    sol_exact = AUX(locations[,1], locations[,2])
    # covariates
    #set.seed(42)
    
    # W=matrix(0,nrow=nobs[i], ncol=1)
    # W[,1]= rnorm(nobs[i], mean=2, sd=0.1)
    
    param =  sol_exact #+ W * beta_ex
    
    mu<-inv.link(param)
    range(mu)
    # sampling response:
    response <- rpois(nobs[i], lambda = mu)
    
    start.GCV = Sys.time()
    output_CPP <- fdaPDE::smooth.FEM(location = locations, 
                                     observations = as.numeric(response), 
                                     FEMbasis =FEMbasis, 
                                     covariates = NULL, #W
                                     max.steps=15, 
                                     fam=FAMILY,
                                     mu0=NULL, 
                                     scale.param=NULL,
                                     lambda = lambda, 
                                     lambda.selection.criterion = 'grid', 
                                     DOF.evaluation = 'exact', 
                                     lambda.selection.lossfunction = 'GCV')
    end.GCV = Sys.time()
    time.GCV = difftime(end.GCV, start.GCV, units="mins")
    print(paste("Time GCV = ", time.GCV,sep=""))
    
    coef.GCV = output_CPP$fit.FEM$coeff[ , output_CPP$optimization$lambda_position]
    mise.GCV = integrate_f( FEM( (coef.ex - coef.GCV)^2, FEMbasis) )
    print(paste("MISE GCV = ", mise.GCV,sep=""))
    mise.GCV
    #R_plot_graph.ggplot2( FEM( coef.GCV, FEMbasis ))
    #R_plot_graph.ggplot2(FEM(coef.ex, FEMbasis))
    
    RMSE.GCV = norm(coef.ex-coef.GCV, type="2")
    print(paste("RMSE GCV = ", RMSE.GCV,sep="")) 
    
    #rmse.GCV.beta = norm(output_CPP$solution$beta[output_CPP$optimization$lambda_position] - beta_ex, type="2")
    #print(paste("RMSE  beta GCV = ", rmse.GCV.beta,sep=""))
  }
}



L = 1125
f = function(x,y){-5/2*sin(2*pi*x/L) * cos(2*pi*y/L) + 4/5*(3*pi*x/L)}
x = seq(0,L, length.out = 40)
y = x
locations = expand.grid(x,y)
sol_exact = f(locations[,1],locations[,2])
contour(locations[,1], locations[,2], sol_exact)
