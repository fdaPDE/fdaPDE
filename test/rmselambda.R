xeval=runif(10000,0,1)
yeval=runif(10000,0,1)
loc_eval=cbind(xeval,yeval)
RMSE<-function(f,g) sqrt(mean((f-g)^2))
rmseD<- NULL
rmseCovD<-NULL

for(l in lambda)
{
  set.seed(1527)
  #best_lambda <- NULL
  rmse <- NULL
  for (i in 1:30)
  {
    output_CPP<-smooth.FEM(observations=data+rnorm (length(nodes[1,]), sd=0.01*(ran[2]-ran[1])),
                           FEMbasis=FEMbasis, #locations=locations, #covariates=cov,
                           lambda=l)
    # best_lambda<-c(best_lambda,
    #               which(output_CPP$optimization$GCV_vector==min(output_CPP$optimization$GCV_vector)))
    # output_CPP<-smooth.FEM(observations=data3d+rnorm (length(nodes3d[1,]), sd=0.01*(ran3d[2]-ran3d[1])),
    #                        FEMbasis=FEMbasis, #locations=locations,
    #                        lambda=lambda[best_lambda[i]])
    rmse <- c(rmse, RMSE(func( loc_eval),eval.FEM(output_CPP$fit.FEM,locations=loc_eval)))
  }
  rmseD<-cbind(rmseD, rmse)
}

for(l in lambda)
{
  set.seed(1527)
  #best_lambda <- NULL
  rmse <- NULL
  for (i in 1:30)
  {
    output_CPP<-smooth.FEM(observations=dataCov+rnorm (length(nodes[1,]), sd=0.01*(ranCov[2]-ranCov[1])),
                           FEMbasis=FEMbasis, covariates=cov,#locations=locations,
                           lambda=l)
    # best_lambda<-c(best_lambda,
    #               which(output_CPP$optimization$GCV_vector==min(output_CPP$optimization$GCV_vector)))
    # output_CPP<-smooth.FEM(observations=data3d+rnorm (length(nodes3d[1,]), sd=0.01*(ran3d[2]-ran3d[1])),
    #                        FEMbasis=FEMbasis, #locations=locations,
    #                        lambda=lambda[best_lambda[i]])
    rmse <- c(rmse, RMSE(func( loc_eval),eval.FEM(output_CPP$fit.FEM,locations=loc_eval)))
  }
  rmseCovD<-cbind(rmseCovD, rmse)
}


rmseGLMA<-NULL
for(l in lambdaGLM)
{
  set.seed(1527)
  #best_lambda <- NULL
  rmse <- NULL
  for (i in 1:30)
  {
    output_CPP<-smooth.FEM(locations = NULL, observations = as.numeric(rpois(nloc, lambda = mu)),
                           FEMbasis =FEMbasis, covariates = desmat,
                           max.steps=15, family =FAMILY, mu0=NULL, scale.param=NULL,
                           lambda = l)
    # best_lambda<-c(best_lambda,
    #               which(output_CPP$optimization$GCV_vector==min(output_CPP$optimization$GCV_vector)))
    # output_CPP<-smooth.FEM(observations=data3d+rnorm (length(nodes3d[1,]), sd=0.01*(ran3d[2]-ran3d[1])),
    #                        FEMbasis=FEMbasis, #locations=locations,
    #                        lambda=lambda[best_lambda[i]])
    rmse <- c(rmse, RMSE(a1* sin(2*pi*loc_eval[,1]) +  a2* sin(2*pi*loc_eval[,2]) +  a3*sin(2*pi*loc_eval[,3]) + 7,eval.FEM(output_CPP$fit.FEM,locations=loc_eval)))
  }
  rmseGLMA<-cbind(rmseGLMA, rmse)
}

set.seed(1527)
#best_lambda <- NULL
for (i in 1:30)
{
  
  
  # best_lambda<-c(best_lambda,
  #               which(output_CPP$optimization$GCV_vector==min(output_CPP$optimization$GCV_vector)))
  # output_CPP<-smooth.FEM(observations=data3d+rnorm (length(nodes3d[1,]), sd=0.01*(ran3d[2]-ran3d[1])),
  #                        FEMbasis=FEMbasis, #locations=locations,
  #                        lambda=lambda[best_lambda[i]])
  rmseA <- c(rmseA, RMSE(fs.test.3D( loc_eval[,1],loc_eval[,2],loc_eval[,3]),eval.FEM(output_CPP$fit.FEM,locations=loc_eval)))
}


set.seed(1527)
#best_lambdaCov <- NULL
for (i in 1:30)
{
  output_CPP<-smooth.FEM(observations=dataCov+rnorm (length(nodes[1,]), sd=0.01*(ranCov[2]-ranCov[1])),
                         FEMbasis=FEMbasis, lambda=lambda, covariates = cov, #locations=locations,
                         lambda.selection.lossfunction = 'GCV')
  best_lambdaCov<-c(best_lambdaCov,
                which(output_CPP$optimization$GCV_vector==min(output_CPP$optimization$GCV_vector)))
  # output_CPP<-smooth.FEM(observations=dataCov+rnorm (length(nodes[1,]), sd=0.01*(ranCov[2]-ranCov[1])),
  #                        FEMbasis=FEMbasis, covariates = cov, #locations=locations,
  #                        lambda=lambda[best_lambdaCov[i]])
  # rmseCovA <- c(rmseCovA, RMSE(func( loc_eval),eval.FEM(output_CPP$fit.FEM,locations=loc_eval)))
}

plot(log10(lambda), colMeans(rmseA), pch=4, col=1, ylab="Mean RMSE")
lines(log10(lambda), colMeans(rmseA), pch=4, col=1)
points(log10(lambda), colMeans(rmseB), pch=4, col=2)
lines(log10(lambda), colMeans(rmseB), col=2)
points(log10(lambda), colMeans(rmseC), pch=1, col=3)
points(log10(lambda), colMeans(rmseD), pch=6, col=4)
grid() 
legend('topleft',legend=c('no preconditioner', 'mass lumping',' lambda preconditioner', 'block preconditioner'), cex=0.62,col=1:4, pch=19)

boxplot(rmseA, names=log10(lambda), ylab="RMSE", xlab="log_10(lambda)")
title("No preconditioner RMSE boxplot")
boxplot(rmseB, names=log10(lambda), ylab="RMSE", xlab="log_10(lambda)")
title("Mass lumping RMSE boxplot")
boxplot(rmseC, names=log10(lambda), ylab="RMSE", xlab="log_10(lambda)")
title("Lambda preconditioner RMSE boxplot")
boxplot(rmseD, names=log10(lambda), ylab="RMSE", xlab="log_10(lambda")
title("Block preconditioner RMSE boxplot")

boxplot(rmseA[,which.min(colMeans(rmseA))],rmseB[,which.min(colMeans(rmseB))],
        rmseC[,which.min(colMeans(rmseC))],rmseD[,which.min(colMeans(rmseD))],
        names=c("fdaPDE", "Mass lumping", "Lambda", "Block"),
        col=c("grey", 2,3,4))

plot(log10(lambda), colMeans(rmseCovA), pch=4, col=1, ylab="Mean RMSE")
lines(log10(lambda), colMeans(rmseCovA), pch=4, col=1)
points(log10(lambda), colMeans(rmseCovB), pch=4, col=2)
lines(log10(lambda), colMeans(rmseCovB), col=2)
points(log10(lambda), colMeans(rmseCovC), pch=1, col=3)
points(log10(lambda), colMeans(rmseCovD), pch=6, col=4)
grid() 
legend('topleft',legend=c('no preconditioner', 'mass lumping',' lambda preconditioner', 'block preconditioner'), cex=0.62,col=1:4, pch=19)

boxplot(rmseCovA, names=log10(lambda), ylab="RMSE", xlab="log_10(lambda)")
title("No preconditioner RMSE boxplot")
boxplot(rmseCovB, names=log10(lambda), ylab="RMSE", xlab="log_10(lambda)")
title("Mass lumping RMSE boxplot")
boxplot(rmseCovC, names=log10(lambda), ylab="RMSE", xlab="log_10(lambda)")
title("Lambda preconditioner RMSE boxplot")
boxplot(rmseCovD, names=log10(lambda), ylab="RMSE", xlab="log_10(lambda)")
title("Block preconditioner RMSE boxplot")

boxplot(rmseCovA[,which.min(colMeans(rmseCovA))],rmseB[,which.min(colMeans(rmseCovB))],
        rmseCovC[,which.min(colMeans(rmseCovC))],rmseD[,which.min(colMeans(rmseCovD))],
        names=c("fdaPDE", "Mass lumping", "Lambda", "Block"),
        col=c("grey", 2,3,4))

plot(log10(lambdaGLM), colMeans(rmseGLMA), pch=4, col=1, ylab="Mean RMSE")
lines(log10(lambdaGLM), colMeans(rmseGLMA), pch=4, col=1)
points(log10(lambdaGLM), colMeans(rmseGLMB), pch=4, col=2)
lines(log10(lambdaGLM), colMeans(rmseGLMB), col=2)
points(log10(lambdaGLM), colMeans(rmseGLMC), pch=4, col=3)
lines(log10(lambdaGLM), colMeans(rmseGLMC), col=3)
points(log10(lambdaGLM), colMeans(rmseGLMD), pch=6, col=4)
lines(log10(lambdaGLM), colMeans(rmseGLMD), col=4)
grid() 
legend('bottomright',legend=c('no preconditioner', 'mass lumping',' lambda preconditioner', 'block preconditioner'), cex=0.62,col=1:4, pch=19)

boxplot(rmseGLMA, names=log10(lambda), ylab="RMSE", xlab="log_10(lambda)")
title("No preconditioner RMSE boxplot")
boxplot(rmseGLMB, names=log10(lambda), ylab="RMSE", xlab="log_10(lambda)")
title("Mass lumping RMSE boxplot")
boxplot(rmseGLMC, names=log10(lambda), ylab="RMSE", xlab="log_10(lambda)")
title("Lambda preconditioner RMSE boxplot")
boxplot(rmseGLMD, names=log10(lambda), ylab="RMSE", xlab="log_10(lambda")
title("Block preconditioner RMSE boxplot")

boxplot(rmseGLMA[,which.min(colMeans(rmseGLMA))],rmseGLMB[,which.min(colMeans(rmseGLMB))],
        rmseGLMC[,which.min(colMeans(rmseGLMC))],rmseGLMD[,which.min(colMeans(rmseGLMD))],
        names=c("fdaPDE", "Mass lumping", "Lambda", "Block"),
        col=c("grey", 2,3,4))


abline(v=log10(mean(best_lambdaCov8[2:20])), lty=3, col=)
