run.R2jags.model <- function(d,
                  nt=10,
                  ni=1100,
                  nb=100,
                  nc=3) {
  
  sink('model.jags')
  cat('model {


	mu~dnorm(0,0.01)
	
	for (sp in 1:nsp) {
		mu.beta.sp[sp]<-mu
	}
	

    lambda.pre ~ dunif(0.0001, 0.9999) #Establish range that lambda can search in
    
    logit.lambda <- logit(lambda.pre)
    logit(lambda) <- logit.lambda
    
    
    beta.mat[1:nsp,1:nsp] <- lambda*VCOV[,] + (1-lambda)*ID[,]
    
    sigma.psi.beta ~ dunif(0,20)
    tau.beta.sp.mat[1:nsp,1:nsp] <- inverse(sigma.psi.beta^2*beta.mat[,])
    ## draw species-specific slopes
    
    
    Y[1:nsp] ~ dmnorm(mu.beta.sp[], tau.beta.sp.mat[,])

  }', fill = TRUE)
  sink()
	
	
	  attach(d$data)
  res <- jags.parallel(data=list('Y', 'nsp', 'VCOV', 'ID'),
                       parameters.to.save=d$params,
                       model.file='model.jags',
                       n.thin=nt,
                       n.iter=ni,
                       n.burnin=nb,
                       n.chains=nc,
                       working.directory=NULL)
  detach(d$data)
  
  res
}

## specify the parameters to be monitored
get.params <- function()
  c('logit.lambda',
  	'lambda',
    'sigma.psi.beta',
    'mu')



	